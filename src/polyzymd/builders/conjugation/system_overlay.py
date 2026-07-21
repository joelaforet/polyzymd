"""OpenMM System merger for mixed GLYCAM/OpenFF conjugate overlays."""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.ownership import OwnershipManifest, TermOwnership

_PARITY_TOLERANCE = 1.0e-10


@dataclass(frozen=True)
class OverlayMergeResult:
    """Result of replacing scoped GLYCAM terms in a baseline OpenMM System."""

    system: Any
    ownership_manifest: OwnershipManifest
    diagnostics: dict[str, Any]

    def save_artifacts(self, output_dir: Path | str) -> dict[str, Path]:
        """Write ownership and overlay diagnostics sidecars."""

        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        ownership_path = self.ownership_manifest.save(directory / "ownership_manifest.json")
        diagnostics_path = directory / "overlay_diagnostics.json"
        diagnostics_path.write_text(
            json.dumps(self.diagnostics, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return {"ownership_manifest": ownership_path, "overlay_diagnostics": diagnostics_path}


def merge_openmm_system_overlay(
    *,
    baseline_system: Any,
    native_system: Any,
    glycam_particles: set[int] | frozenset[int],
    atom_mapping: dict[int, int] | None = None,
    attachments: tuple[dict[str, Any], ...] = (),
) -> OverlayMergeResult:
    """Replace GLYCAM-owned terms in a complete baseline System.

    ``atom_mapping`` maps native reference atom indices to final baseline atom
    indices. It may be a strict subset so a GLYCAM-only native reference can
    overlay a full OpenFF baseline that also contains Sage polymer, solvent, and
    ions.
    """

    baseline_count = baseline_system.getNumParticles()
    native_count = native_system.getNumParticles()
    mapping = atom_mapping or _identity_mapping_for_equal_counts(baseline_count, native_count)
    _validate_subset_mapping(mapping, baseline_count, native_count)
    native_to_baseline = dict(mapping)
    baseline_to_native = {baseline: native for native, baseline in native_to_baseline.items()}
    native_glycam = set(glycam_particles)
    _validate_particle_indices(native_glycam, native_count, label="native GLYCAM")
    unmapped_glycam = sorted(native_glycam - set(native_to_baseline))
    if unmapped_glycam:
        raise ValueError(
            f"Native GLYCAM particles are missing baseline mappings: {unmapped_glycam}"
        )
    glycam = {native_to_baseline[index] for index in native_glycam}

    system = _empty_system_like(baseline_system)
    terms: list[TermOwnership] = []
    diagnostics: dict[str, Any] = {
        "route": "mixed_overlay",
        "mapping_count": len(native_to_baseline),
        "glycam_baseline_indices": sorted(glycam),
        "forces": {},
        "conflicts": [],
    }
    _copy_particles(system, baseline_system, native_system, glycam, baseline_to_native, terms)
    _copy_constraints(
        system,
        baseline_system,
        native_system,
        glycam,
        native_glycam,
        native_to_baseline,
        terms,
        diagnostics,
    )
    _copy_forces(
        system,
        baseline_system,
        native_system,
        glycam,
        native_glycam,
        native_to_baseline,
        baseline_to_native,
        terms,
        diagnostics,
    )
    diagnostics["counts"] = _system_counts(system)
    diagnostics["parity"] = _audit_overlay_parity(
        merged_system=system,
        baseline_system=baseline_system,
        native_system=native_system,
        glycam=glycam,
        native_glycam=native_glycam,
        native_to_baseline=native_to_baseline,
        baseline_to_native=baseline_to_native,
    )
    manifest = OwnershipManifest(
        domains={
            "glycam": sorted(glycam),
            "generic": [index for index in range(baseline_count) if index not in glycam],
        },
        attachments=list(attachments),
        atom_mapping=native_to_baseline,
        terms=tuple(terms),
        conflicts=tuple(diagnostics["conflicts"]),
    )
    manifest.validate_complete(baseline_count)
    return OverlayMergeResult(system=system, ownership_manifest=manifest, diagnostics=diagnostics)


def _empty_system_like(source: Any) -> Any:
    """Create an empty OpenMM System with source-level box settings."""

    import openmm

    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(*source.getDefaultPeriodicBoxVectors())
    return system


def _copy_particles(
    target: Any,
    baseline: Any,
    native: Any,
    glycam: set[int],
    baseline_to_native: dict[int, int],
    terms: list[TermOwnership],
) -> None:
    """Copy particles, choosing native masses for mapped GLYCAM particles."""

    for index in range(baseline.getNumParticles()):
        owner = "glycam" if index in glycam else "generic"
        if owner == "glycam":
            native_index = baseline_to_native[index]
            target.addParticle(native.getParticleMass(native_index))
            if native.isVirtualSite(native_index):
                raise ValueError("Overlay does not support remapped native virtual sites")
        else:
            target.addParticle(baseline.getParticleMass(index))
            if baseline.isVirtualSite(index):
                target.setVirtualSite(index, baseline.getVirtualSite(index))
        terms.append(TermOwnership(kind="particle", owner=owner, atoms=(index,), source=owner))


def _copy_constraints(
    target: Any,
    baseline: Any,
    native: Any,
    glycam: set[int],
    native_glycam: set[int],
    native_to_baseline: dict[int, int],
    terms: list[TermOwnership],
    diagnostics: dict[str, Any],
) -> None:
    """Copy constraints by ownership."""

    kept = _add_constraints_from(target, baseline, glycam, False, terms, "generic")
    skipped = _count_constraints_touching(baseline, glycam)
    added = _add_constraints_from(
        target,
        native,
        native_glycam,
        True,
        terms,
        "glycam",
        atom_mapping=native_to_baseline,
    )
    diagnostics["constraints"] = {
        "kept_baseline": kept,
        "skipped_baseline_touching_glycam": skipped,
        "added_native_glycam": added,
    }


def _copy_forces(
    target: Any,
    baseline: Any,
    native: Any,
    glycam: set[int],
    native_glycam: set[int],
    native_to_baseline: dict[int, int],
    baseline_to_native: dict[int, int],
    terms: list[TermOwnership],
    diagnostics: dict[str, Any],
) -> None:
    """Copy supported forces with scoped overlay and preserve unrelated forces."""

    handled = ("NonbondedForce", "HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce")
    baseline_forces = _forces_by_name(baseline)
    native_forces = _forces_by_name(native)
    for name in handled:
        if name not in baseline_forces or name not in native_forces:
            continue
        if name == "NonbondedForce":
            force, force_terms, force_diag = _merged_nonbonded_force(
                baseline_forces[name],
                native_forces[name],
                glycam,
                native_glycam,
                native_to_baseline,
                baseline_to_native,
            )
        else:
            kind = {"HarmonicBondForce": "bond", "HarmonicAngleForce": "angle"}.get(name, "torsion")
            force, force_terms, force_diag = _merged_indexed_force(
                baseline_forces[name],
                native_forces[name],
                glycam,
                native_glycam,
                native_to_baseline,
                kind,
            )
        _copy_force_metadata(force, baseline_forces[name])
        force_diag["force_metadata"] = _force_metadata_audit(force)
        target.addForce(force)
        terms.extend(force_terms)
        diagnostics["forces"][name] = force_diag

    for force in baseline.getForces():
        class_name = force.__class__.__name__
        if class_name not in handled:
            policy = _baseline_unhandled_force_policy(force, glycam)
            diagnostics["forces"][class_name] = policy
            if policy["action"] == "reject":
                raise ValueError(policy["message"])
            target.addForce(copy.deepcopy(force))


def _merged_nonbonded_force(
    baseline_force: Any,
    native_force: Any,
    glycam: set[int],
    native_glycam: set[int],
    native_to_baseline: dict[int, int],
    baseline_to_native: dict[int, int],
) -> tuple[Any, list[TermOwnership], dict[str, int]]:
    """Merge NonbondedForce particles and exceptions."""

    import openmm

    force = openmm.NonbondedForce()
    _copy_nonbonded_settings(force, baseline_force)
    terms: list[TermOwnership] = []
    for index in range(baseline_force.getNumParticles()):
        if index in glycam:
            force.addParticle(*native_force.getParticleParameters(baseline_to_native[index]))
        else:
            force.addParticle(*baseline_force.getParticleParameters(index))
    kept = _add_exceptions_from(force, baseline_force, glycam, False, terms, "generic")
    skipped = _count_exceptions_touching(baseline_force, glycam)
    added = _add_exceptions_from(
        force,
        native_force,
        native_glycam,
        True,
        terms,
        "glycam",
        atom_mapping=native_to_baseline,
    )
    return (
        force,
        terms,
        {
            "kept_baseline_exceptions": kept,
            "skipped_baseline_touching_glycam": skipped,
            "added_native_glycam_exceptions": added,
            "global_settings": _nonbonded_settings_audit(force),
            "native_global_settings": _nonbonded_settings_audit(native_force),
        },
    )


def _merged_indexed_force(
    baseline_force: Any,
    native_force: Any,
    glycam: set[int],
    native_glycam: set[int],
    native_to_baseline: dict[int, int],
    kind: str,
) -> tuple[Any, list[TermOwnership], dict[str, int]]:
    """Merge bond-like forces using their indexed term getters."""

    force = _new_force_like(baseline_force)
    terms: list[TermOwnership] = []
    kept = _add_terms_from(force, baseline_force, glycam, False, terms, "generic", kind)
    skipped = _count_terms_touching(baseline_force, glycam, kind)
    added = _add_terms_from(
        force,
        native_force,
        native_glycam,
        True,
        terms,
        "glycam",
        kind,
        atom_mapping=native_to_baseline,
    )
    return (
        force,
        terms,
        {
            "kept_baseline": kept,
            "skipped_baseline_touching_glycam": skipped,
            "added_native_glycam": added,
        },
    )


def _new_force_like(source: Any) -> Any:
    """Create an empty force of the same OpenMM class."""

    import openmm

    force = getattr(openmm, source.__class__.__name__)()
    _copy_force_metadata(force, source)
    return force


def _copy_force_metadata(target: Any, source: Any) -> None:
    """Copy force-level metadata available through the OpenMM API."""

    if hasattr(target, "setName") and hasattr(source, "getName"):
        target.setName(source.getName())
    if hasattr(target, "setForceGroup") and hasattr(source, "getForceGroup"):
        target.setForceGroup(source.getForceGroup())
    if hasattr(target, "setUsesPeriodicBoundaryConditions") and hasattr(
        source, "usesPeriodicBoundaryConditions"
    ):
        target.setUsesPeriodicBoundaryConditions(source.usesPeriodicBoundaryConditions())


def _force_metadata_audit(force: Any) -> dict[str, Any]:
    """Return JSON-safe force-level metadata."""

    return {
        "name": force.getName() if hasattr(force, "getName") else "",
        "force_group": force.getForceGroup() if hasattr(force, "getForceGroup") else None,
        "uses_pbc": (
            force.usesPeriodicBoundaryConditions()
            if hasattr(force, "usesPeriodicBoundaryConditions")
            else None
        ),
    }


def _add_constraints_from(
    target: Any,
    source: Any,
    glycam: set[int],
    include_touching: bool,
    terms: list[TermOwnership],
    owner: str,
    atom_mapping: dict[int, int] | None = None,
) -> int:
    """Add constraints matching the ownership predicate."""

    count = 0
    for index in range(source.getNumConstraints()):
        i, j, distance = source.getConstraintParameters(index)
        touching = i in glycam or j in glycam
        if touching is include_touching:
            atoms = _map_atoms((i, j), atom_mapping)
            target.addConstraint(*atoms, distance)
            terms.append(TermOwnership(kind="constraint", owner=owner, atoms=atoms, source=owner))
            count += 1
    return count


def _add_exceptions_from(
    target: Any,
    source: Any,
    glycam: set[int],
    include_touching: bool,
    terms: list[TermOwnership],
    owner: str,
    atom_mapping: dict[int, int] | None = None,
) -> int:
    """Add exceptions matching the ownership predicate."""

    count = 0
    for index in range(source.getNumExceptions()):
        i, j, chargeprod, sigma, epsilon = source.getExceptionParameters(index)
        touching = i in glycam or j in glycam
        if touching is include_touching:
            atoms = _map_atoms((i, j), atom_mapping)
            target.addException(*atoms, chargeprod, sigma, epsilon, replace=False)
            terms.append(TermOwnership(kind="exception", owner=owner, atoms=atoms, source=owner))
            count += 1
    return count


def _add_terms_from(
    target: Any,
    source: Any,
    glycam: set[int],
    include_touching: bool,
    terms: list[TermOwnership],
    owner: str,
    kind: str,
    atom_mapping: dict[int, int] | None = None,
) -> int:
    """Add bond, angle, or torsion terms matching the ownership predicate."""

    count = 0
    getter_name = {
        "bond": "getBondParameters",
        "angle": "getAngleParameters",
        "torsion": "getTorsionParameters",
    }[kind]
    adder_name = {"bond": "addBond", "angle": "addAngle", "torsion": "addTorsion"}[kind]
    count_name = {"bond": "getNumBonds", "angle": "getNumAngles", "torsion": "getNumTorsions"}[kind]
    atom_count = {"bond": 2, "angle": 3, "torsion": 4}[kind]
    for index in range(getattr(source, count_name)()):
        params = getattr(source, getter_name)(index)
        atoms = tuple(int(value) for value in params[:atom_count])
        touching = any(atom in glycam for atom in atoms)
        if touching is include_touching:
            mapped_atoms = _map_atoms(atoms, atom_mapping)
            getattr(target, adder_name)(*mapped_atoms, *params[atom_count:])
            terms.append(TermOwnership(kind=kind, owner=owner, atoms=mapped_atoms, source=owner))
            count += 1
    return count


def _copy_nonbonded_settings(target: Any, source: Any) -> None:
    """Copy global NonbondedForce settings."""

    target.setNonbondedMethod(source.getNonbondedMethod())
    target.setCutoffDistance(source.getCutoffDistance())
    target.setUseSwitchingFunction(source.getUseSwitchingFunction())
    target.setSwitchingDistance(source.getSwitchingDistance())
    target.setUseDispersionCorrection(source.getUseDispersionCorrection())
    target.setEwaldErrorTolerance(source.getEwaldErrorTolerance())
    if hasattr(target, "setPMEParameters") and hasattr(source, "getPMEParameters"):
        target.setPMEParameters(*source.getPMEParameters())
    if hasattr(target, "setLJPMEParameters") and hasattr(source, "getLJPMEParameters"):
        target.setLJPMEParameters(*source.getLJPMEParameters())
    if hasattr(target, "setReciprocalSpaceForceGroup") and hasattr(
        source, "getReciprocalSpaceForceGroup"
    ):
        target.setReciprocalSpaceForceGroup(source.getReciprocalSpaceForceGroup())
    _copy_force_metadata(target, source)


def _nonbonded_settings_audit(force: Any) -> dict[str, Any]:
    """Return JSON-safe global NonbondedForce settings."""

    audit = {
        **_force_metadata_audit(force),
        "method": int(force.getNonbondedMethod()),
        "cutoff_distance": _canonical_value(force.getCutoffDistance()),
        "use_switching_function": force.getUseSwitchingFunction(),
        "switching_distance": _canonical_value(force.getSwitchingDistance()),
        "use_dispersion_correction": force.getUseDispersionCorrection(),
        "ewald_error_tolerance": float(force.getEwaldErrorTolerance()),
    }
    if hasattr(force, "getPMEParameters"):
        alpha, nx, ny, nz = force.getPMEParameters()
        audit["pme_parameters"] = [_canonical_value(alpha), int(nx), int(ny), int(nz)]
    if hasattr(force, "getLJPMEParameters"):
        alpha, nx, ny, nz = force.getLJPMEParameters()
        audit["ljpme_parameters"] = [_canonical_value(alpha), int(nx), int(ny), int(nz)]
    if hasattr(force, "getReciprocalSpaceForceGroup"):
        audit["reciprocal_space_force_group"] = int(force.getReciprocalSpaceForceGroup())
    return audit


def _baseline_unhandled_force_policy(force: Any, glycam: set[int]) -> dict[str, Any]:
    """Return whether an unhandled baseline force can be copied safely."""

    class_name = force.__class__.__name__
    touched = _force_touched_particles(force)
    if touched is None:
        if class_name.startswith("Custom"):
            return {
                "action": "reject",
                "reason": "non_introspectable_custom_force",
                "message": (
                    "Mixed overlay refuses to copy baseline custom force "
                    f"{class_name!r} because its particle ownership cannot be audited"
                ),
            }
        return {"action": "copy", "reason": "non_custom_force", "touched_particle_count": None}
    overlap = sorted(touched & glycam)
    if overlap:
        return {
            "action": "reject",
            "reason": "custom_force_touches_glycam",
            "touched_particle_count": len(touched),
            "glycam_overlap": overlap,
            "message": (
                "Mixed overlay refuses to copy baseline custom force "
                f"{class_name!r} because it touches GLYCAM-owned atoms {overlap[:10]}"
            ),
        }
    return {
        "action": "copy",
        "reason": "generic_only_or_global_force",
        "touched_particle_count": len(touched),
    }


def _force_touched_particles(force: Any) -> set[int] | None:
    """Return particles referenced by a force, or ``None`` when unknown."""

    class_name = force.__class__.__name__
    if class_name == "CustomExternalForce":
        return _atoms_from_external_force(force)
    if class_name == "CustomCompoundBondForce":
        return _atoms_from_compound_bond_force(force)
    if class_name == "CustomCentroidBondForce":
        return _atoms_from_centroid_force(force)
    if class_name == "CustomCVForce":
        return _atoms_from_cv_force(force)
    if hasattr(force, "getNumParticles") and hasattr(force, "getParticleParameters"):
        return set(range(force.getNumParticles()))
    if hasattr(force, "getNumBonds") and hasattr(force, "getBondParameters"):
        return _atoms_from_indexed_parameters(force, "getNumBonds", "getBondParameters", 2)
    if hasattr(force, "getNumAngles") and hasattr(force, "getAngleParameters"):
        return _atoms_from_indexed_parameters(force, "getNumAngles", "getAngleParameters", 3)
    if hasattr(force, "getNumTorsions") and hasattr(force, "getTorsionParameters"):
        return _atoms_from_indexed_parameters(force, "getNumTorsions", "getTorsionParameters", 4)
    return set() if not class_name.startswith("Custom") else None


def _atoms_from_indexed_parameters(
    force: Any, count_name: str, getter_name: str, atom_count: int
) -> set[int]:
    """Return atom indices from fixed-width indexed force parameters."""

    atoms: set[int] = set()
    for index in range(getattr(force, count_name)()):
        params = getattr(force, getter_name)(index)
        atoms.update(int(value) for value in params[:atom_count])
    return atoms


def _atoms_from_external_force(force: Any) -> set[int]:
    """Return atom indices referenced by a CustomExternalForce."""

    return {int(force.getParticleParameters(index)[0]) for index in range(force.getNumParticles())}


def _atoms_from_compound_bond_force(force: Any) -> set[int]:
    """Return atom indices referenced by a CustomCompoundBondForce."""

    atoms: set[int] = set()
    for index in range(force.getNumBonds()):
        particles = force.getBondParameters(index)[0]
        atoms.update(int(value) for value in particles)
    return atoms


def _atoms_from_centroid_force(force: Any) -> set[int]:
    """Return atom indices referenced by a CustomCentroidBondForce."""

    atoms: set[int] = set()
    for index in range(force.getNumGroups()):
        particles = force.getGroupParameters(index)[0]
        atoms.update(int(value) for value in particles)
    return atoms


def _atoms_from_cv_force(force: Any) -> set[int] | None:
    """Return atom indices referenced by nested CustomCVForce variables."""

    atoms: set[int] = set()
    for index in range(force.getNumCollectiveVariables()):
        nested = force.getCollectiveVariable(index)
        nested_atoms = _force_touched_particles(nested)
        if nested_atoms is None:
            return None
        atoms.update(nested_atoms)
    return atoms


def _forces_by_name(system: Any) -> dict[str, Any]:
    """Return unique forces by class name or fail on duplicates."""

    forces: dict[str, Any] = {}
    for force in system.getForces():
        name = force.__class__.__name__
        if name in forces:
            raise ValueError(f"Overlay does not support multiple {name} instances")
        forces[name] = force
    return forces


def _map_atoms(atoms: tuple[int, ...], atom_mapping: dict[int, int] | None) -> tuple[int, ...]:
    """Map native atom indices to baseline atom indices when requested."""

    if atom_mapping is None:
        return atoms
    missing = [atom for atom in atoms if atom not in atom_mapping]
    if missing:
        raise ValueError(f"Native GLYCAM term contains unmapped atoms: {missing}")
    return tuple(atom_mapping[atom] for atom in atoms)


def _identity_mapping_for_equal_counts(
    baseline_particle_count: int, native_particle_count: int
) -> dict[int, int]:
    """Return identity mapping only when Systems have identical particle counts."""

    if baseline_particle_count != native_particle_count:
        raise ValueError("Overlay atom mapping is required when particle counts differ")
    return {index: index for index in range(baseline_particle_count)}


def _validate_subset_mapping(
    mapping: dict[int, int], baseline_particle_count: int, native_particle_count: int
) -> None:
    """Validate native-to-baseline subset mapping for an overlay reference."""

    _validate_particle_indices(set(mapping), native_particle_count, label="native mapped")
    _validate_particle_indices(
        set(mapping.values()), baseline_particle_count, label="baseline mapped"
    )
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(
            "Overlay atom mapping must not map multiple native atoms to one baseline atom"
        )


def _count_constraints_touching(source: Any, domain: set[int]) -> int:
    """Count constraints touching a domain."""

    return sum(
        1
        for index in range(source.getNumConstraints())
        if any(atom in domain for atom in source.getConstraintParameters(index)[:2])
    )


def _count_exceptions_touching(source: Any, domain: set[int]) -> int:
    """Count exceptions touching a domain."""

    return sum(
        1
        for index in range(source.getNumExceptions())
        if any(atom in domain for atom in source.getExceptionParameters(index)[:2])
    )


def _count_terms_touching(source: Any, domain: set[int], kind: str) -> int:
    """Count bond-like terms touching a domain."""

    getter_name = {
        "bond": "getBondParameters",
        "angle": "getAngleParameters",
        "torsion": "getTorsionParameters",
    }[kind]
    count_name = {"bond": "getNumBonds", "angle": "getNumAngles", "torsion": "getNumTorsions"}[kind]
    atom_count = {"bond": 2, "angle": 3, "torsion": 4}[kind]
    return sum(
        1
        for index in range(getattr(source, count_name)())
        if any(atom in domain for atom in getattr(source, getter_name)(index)[:atom_count])
    )


def _system_counts(system: Any) -> dict[str, int]:
    """Return JSON-safe force and particle counts for diagnostics."""

    counts = {"particles": system.getNumParticles(), "constraints": system.getNumConstraints()}
    for force in system.getForces():
        name = force.__class__.__name__
        if name == "NonbondedForce":
            counts["nb_particles"] = force.getNumParticles()
            counts["nb_exceptions"] = force.getNumExceptions()
            counts["nb_method"] = int(force.getNonbondedMethod())
        elif name == "HarmonicBondForce":
            counts["bonds"] = force.getNumBonds()
        elif name == "HarmonicAngleForce":
            counts["angles"] = force.getNumAngles()
        elif name == "PeriodicTorsionForce":
            counts["torsions"] = force.getNumTorsions()
        else:
            counts[name] = counts.get(name, 0) + 1
    return counts


def _audit_overlay_parity(
    *,
    merged_system: Any,
    baseline_system: Any,
    native_system: Any,
    glycam: set[int],
    native_glycam: set[int],
    native_to_baseline: dict[int, int],
    baseline_to_native: dict[int, int],
) -> dict[str, Any]:
    """Validate merged parameters against their selected source Systems."""

    baseline_forces = _forces_by_name(baseline_system)
    native_forces = _forces_by_name(native_system)
    merged_forces = _forces_by_name(merged_system)
    generic = {index for index in range(baseline_system.getNumParticles()) if index not in glycam}
    parity = {
        "tolerance": _PARITY_TOLERANCE,
        "unmapped_glycam_count": len(native_glycam - set(native_to_baseline)),
        "unmapped_generic_count": len(glycam - set(baseline_to_native)),
        "conflict_count": 0,
        "glycam_vs_native": {},
        "generic_vs_baseline": {},
    }
    parity["glycam_vs_native"]["particles"] = _audit_particles(
        merged_system,
        native_system,
        merged_forces["NonbondedForce"],
        native_forces["NonbondedForce"],
        sorted(glycam),
        {baseline: native for native, baseline in native_to_baseline.items()},
        "GLYCAM particle",
    )
    parity["generic_vs_baseline"]["particles"] = _audit_particles(
        merged_system,
        baseline_system,
        merged_forces["NonbondedForce"],
        baseline_forces["NonbondedForce"],
        sorted(generic),
        {index: index for index in generic},
        "generic particle",
    )
    parity["glycam_vs_native"]["constraints"] = _audit_constraint_terms(
        merged_system,
        native_system,
        native_glycam,
        include_touching=True,
        atom_mapping=native_to_baseline,
        label="GLYCAM constraints",
    )
    parity["generic_vs_baseline"]["constraints"] = _audit_constraint_terms(
        merged_system,
        baseline_system,
        glycam,
        include_touching=False,
        atom_mapping=None,
        label="generic constraints",
    )
    for force_name, kind in (
        ("HarmonicBondForce", "bond"),
        ("HarmonicAngleForce", "angle"),
        ("PeriodicTorsionForce", "torsion"),
    ):
        parity["glycam_vs_native"][kind] = _audit_indexed_terms(
            merged_forces[force_name],
            native_forces[force_name],
            native_glycam,
            include_touching=True,
            atom_mapping=native_to_baseline,
            kind=kind,
            label=f"GLYCAM {kind} terms",
        )
        parity["generic_vs_baseline"][kind] = _audit_indexed_terms(
            merged_forces[force_name],
            baseline_forces[force_name],
            glycam,
            include_touching=False,
            atom_mapping=None,
            kind=kind,
            label=f"generic {kind} terms",
        )
    parity["glycam_vs_native"]["exceptions"] = _audit_exception_terms(
        merged_forces["NonbondedForce"],
        native_forces["NonbondedForce"],
        native_glycam,
        include_touching=True,
        atom_mapping=native_to_baseline,
        label="GLYCAM exceptions",
    )
    parity["generic_vs_baseline"]["exceptions"] = _audit_exception_terms(
        merged_forces["NonbondedForce"],
        baseline_forces["NonbondedForce"],
        glycam,
        include_touching=False,
        atom_mapping=None,
        label="generic exceptions",
    )
    return parity


def _audit_particles(
    merged_system: Any,
    source_system: Any,
    merged_nonbonded: Any,
    source_nonbonded: Any,
    merged_indices: list[int],
    merged_to_source: dict[int, int],
    label: str,
) -> dict[str, Any]:
    """Audit particle mass and nonbonded parameter parity."""

    max_delta = 0.0
    components = {"mass": 0.0, "charge": 0.0, "sigma": 0.0, "epsilon": 0.0}
    for merged_index in merged_indices:
        source_index = merged_to_source.get(merged_index)
        if source_index is None:
            raise ValueError(
                f"{label} parity failed: missing source mapping for atom {merged_index}"
            )
        merged_values = (
            _canonical_value(merged_system.getParticleMass(merged_index)),
            *(
                _canonical_value(value)
                for value in merged_nonbonded.getParticleParameters(merged_index)
            ),
        )
        source_values = (
            _canonical_value(source_system.getParticleMass(source_index)),
            *(
                _canonical_value(value)
                for value in source_nonbonded.getParticleParameters(source_index)
            ),
        )
        for name, merged_value, source_value in zip(
            components, merged_values, source_values, strict=True
        ):
            delta = abs(merged_value - source_value)
            components[name] = max(components[name], delta)
            max_delta = max(max_delta, delta)
    _fail_on_delta(label, max_delta)
    return {
        "count": len(merged_indices),
        "missing_count": 0,
        "extra_count": 0,
        "max_abs_delta": max_delta,
        "component_max_abs_deltas": components,
    }


def _audit_constraint_terms(
    merged_system: Any,
    source_system: Any,
    domain: set[int],
    *,
    include_touching: bool,
    atom_mapping: dict[int, int] | None,
    label: str,
) -> dict[str, Any]:
    """Audit constraint multiset parity."""

    merged = _constraint_records(
        merged_system, _map_domain(domain, atom_mapping), include_touching, None
    )
    expected = _constraint_records(source_system, domain, include_touching, atom_mapping)
    return _compare_record_multisets(label, merged, expected)


def _audit_exception_terms(
    merged_force: Any,
    source_force: Any,
    domain: set[int],
    *,
    include_touching: bool,
    atom_mapping: dict[int, int] | None,
    label: str,
) -> dict[str, Any]:
    """Audit NonbondedForce exception multiset parity."""

    merged = _exception_records(
        merged_force, _map_domain(domain, atom_mapping), include_touching, None
    )
    expected = _exception_records(source_force, domain, include_touching, atom_mapping)
    return _compare_record_multisets(label, merged, expected)


def _audit_indexed_terms(
    merged_force: Any,
    source_force: Any,
    domain: set[int],
    *,
    include_touching: bool,
    atom_mapping: dict[int, int] | None,
    kind: str,
    label: str,
) -> dict[str, Any]:
    """Audit bond-like force multiset parity."""

    merged = _indexed_records(
        merged_force, _map_domain(domain, atom_mapping), include_touching, None, kind
    )
    expected = _indexed_records(source_force, domain, include_touching, atom_mapping, kind)
    return _compare_record_multisets(label, merged, expected)


def _constraint_records(
    system: Any,
    domain: set[int],
    include_touching: bool,
    atom_mapping: dict[int, int] | None,
) -> dict[tuple[int, ...], list[tuple[float, ...]]]:
    """Return selected constraint records keyed by canonical atom tuple."""

    records: dict[tuple[int, ...], list[tuple[float, ...]]] = {}
    for index in range(system.getNumConstraints()):
        i, j, distance = system.getConstraintParameters(index)
        if (i in domain or j in domain) is include_touching:
            atoms = _canonical_atom_key(_map_atoms((i, j), atom_mapping), "constraint")
            records.setdefault(atoms, []).append((_canonical_value(distance),))
    return records


def _exception_records(
    force: Any,
    domain: set[int],
    include_touching: bool,
    atom_mapping: dict[int, int] | None,
) -> dict[tuple[int, ...], list[tuple[float, ...]]]:
    """Return selected exception records keyed by canonical atom tuple."""

    records: dict[tuple[int, ...], list[tuple[float, ...]]] = {}
    for index in range(force.getNumExceptions()):
        i, j, chargeprod, sigma, epsilon = force.getExceptionParameters(index)
        if (i in domain or j in domain) is include_touching:
            atoms = _canonical_atom_key(_map_atoms((i, j), atom_mapping), "exception")
            records.setdefault(atoms, []).append(
                tuple(_canonical_value(value) for value in (chargeprod, sigma, epsilon))
            )
    return records


def _indexed_records(
    force: Any,
    domain: set[int],
    include_touching: bool,
    atom_mapping: dict[int, int] | None,
    kind: str,
) -> dict[tuple[int, ...], list[tuple[float, ...]]]:
    """Return selected bond-like records keyed by canonical atom tuple."""

    getter_name = {
        "bond": "getBondParameters",
        "angle": "getAngleParameters",
        "torsion": "getTorsionParameters",
    }[kind]
    count_name = {"bond": "getNumBonds", "angle": "getNumAngles", "torsion": "getNumTorsions"}[kind]
    atom_count = {"bond": 2, "angle": 3, "torsion": 4}[kind]
    records: dict[tuple[int, ...], list[tuple[float, ...]]] = {}
    for index in range(getattr(force, count_name)()):
        params = getattr(force, getter_name)(index)
        atoms = tuple(int(value) for value in params[:atom_count])
        if any(atom in domain for atom in atoms) is include_touching:
            key = _canonical_atom_key(_map_atoms(atoms, atom_mapping), kind)
            records.setdefault(key, []).append(
                tuple(_canonical_value(value) for value in params[atom_count:])
            )
    return records


def _compare_record_multisets(
    label: str,
    observed: dict[tuple[int, ...], list[tuple[float, ...]]],
    expected: dict[tuple[int, ...], list[tuple[float, ...]]],
) -> dict[str, Any]:
    """Compare selected parameter records and fail on mismatches."""

    observed_keys = set(observed)
    expected_keys = set(expected)
    missing = sorted(expected_keys - observed_keys)
    extra = sorted(observed_keys - expected_keys)
    if missing or extra:
        raise ValueError(f"{label} parity failed: missing={missing[:5]}, extra={extra[:5]}")
    max_delta = 0.0
    count = 0
    for key in sorted(expected):
        observed_values = sorted(observed[key])
        expected_values = sorted(expected[key])
        if len(observed_values) != len(expected_values):
            raise ValueError(
                f"{label} parity failed for atoms {key}: observed {len(observed_values)} terms, "
                f"expected {len(expected_values)} terms"
            )
        count += len(expected_values)
        for observed_record, expected_record in zip(observed_values, expected_values, strict=True):
            if len(observed_record) != len(expected_record):
                raise ValueError(f"{label} parity failed for atoms {key}: parameter width mismatch")
            record_delta = max(
                (
                    abs(left - right)
                    for left, right in zip(observed_record, expected_record, strict=True)
                ),
                default=0.0,
            )
            max_delta = max(max_delta, record_delta)
    _fail_on_delta(label, max_delta)
    return {
        "count": count,
        "missing_count": 0,
        "extra_count": 0,
        "max_abs_delta": max_delta,
    }


def _canonical_atom_key(atoms: tuple[int, ...], kind: str) -> tuple[int, ...]:
    """Return an atom tuple key with reversible OpenMM terms normalized."""

    if kind in {"bond", "constraint", "exception"}:
        return tuple(sorted(atoms))
    reversed_atoms = tuple(reversed(atoms))
    return min(atoms, reversed_atoms)


def _map_domain(domain: set[int], atom_mapping: dict[int, int] | None) -> set[int]:
    """Map a source domain into merged atom indices when needed."""

    if atom_mapping is None:
        return set(domain)
    return {atom_mapping[index] for index in domain if index in atom_mapping}


def _canonical_value(value: Any) -> float:
    """Return a numeric value in OpenMM's MD unit system when applicable."""

    try:
        from openmm import unit

        return float(value.value_in_unit_system(unit.md_unit_system))
    except AttributeError:
        return float(value)


def _fail_on_delta(label: str, max_delta: float) -> None:
    """Raise when a parity delta exceeds the tight overlay tolerance."""

    if max_delta > _PARITY_TOLERANCE:
        raise ValueError(
            f"{label} parity failed: max absolute parameter delta {max_delta:.3e} exceeds "
            f"{_PARITY_TOLERANCE:.3e}"
        )


def _validate_particle_indices(
    indices: set[int], particle_count: int, *, label: str = "overlay"
) -> None:
    """Validate zero-based particle indices."""

    invalid = sorted(index for index in indices if index < 0 or index >= particle_count)
    if invalid:
        raise ValueError(f"{label} particle indices are out of range: {invalid}")
