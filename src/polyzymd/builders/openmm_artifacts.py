"""Helpers for building and validating OpenMM handoff artifacts."""

from __future__ import annotations

import filecmp
import hashlib
import json
import math
import shutil
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.structure.inspection import (
    COMMON_ION_RESIDUES,
    COMMON_SOLVENT_RESIDUES,
    WATER_RESIDUES,
)
from polyzymd.builders.conjugation.structure.parsing import parse_pdb_conect_adjacency
from polyzymd.config.schema import BuildScope

SOLVATED_SYSTEM_PDB = "solvated_system.pdb"
SOLVATED_CONJUGATE_PDB = "solvated_conjugate_free_polymers.pdb"
_RESTRAINT_FORCE_PREFIX = "PolyzyMD restraint:"
_RESTRAINT_CUSTOM_BOND_EXPRESSIONS = {
    "step(r - r0) * 0.5 * k * (r - r0)^2",
    "step(r0 - r) * 0.5 * k * (r0 - r)^2",
}


@dataclass(frozen=True)
class OpenMMBuildArtifacts:
    """Container describing a completed OpenMM build handoff.

    The container preserves route-specific builder and result objects for
    export callers while standardizing the files required by OpenMM simulation
    startup.
    """

    builder: Any
    interchange: Any | None
    result: Any | None
    pdb_path: Path
    system_path: Path | None
    conjugation_enabled: bool
    scope: BuildScope = BuildScope.SYSTEM
    exact_export_bundle: Any | None = None

    def get_component_info(self) -> Any:
        """Return component metadata from the route-specific build object.

        Returns
        -------
        Any
            Component metadata suitable for exporter calls.
        """
        if self.result is not None:
            if self.scope is BuildScope.STRUCTURE:
                raise RuntimeError("Structure-scope builds do not contain component parameters.")
            return self.result.get_component_info()
        return self.builder.get_component_info()

    def require_final_interchange(self) -> Any:
        """Return the final Interchange from the route-specific build object.

        Returns
        -------
        Any
            Final Interchange object for export.
        """
        if self.scope is BuildScope.STRUCTURE:
            raise RuntimeError(
                "Structure-scope builds do not contain a final Interchange. "
                "Use build.scope: solute or system before exporting."
            )
        if self.result is not None:
            return self.result.require_final_interchange()
        if self.exact_export_bundle is not None:
            return self.exact_export_bundle
        return self.interchange


def conjugation_enabled(sim_config: object) -> bool:
    """Return whether a config requests the conjugation build workflow.

    Parameters
    ----------
    sim_config : object
        Simulation configuration with an optional ``conjugation`` attribute.

    Returns
    -------
    bool
        ``True`` when conjugation is configured and enabled.
    """
    conjugation = getattr(sim_config, "conjugation", None)
    return conjugation is not None and getattr(conjugation, "enabled", False) is True


def ensure_conjugation_pdb_alias(working_dir: Path) -> Path:
    """Create or refresh the standard OpenMM PDB alias for conjugation builds.

    Parameters
    ----------
    working_dir : Path
        Directory containing conjugation build outputs.

    Returns
    -------
    Path
        Path to ``solvated_system.pdb``.

    Raises
    ------
    FileNotFoundError
        If the route-specific conjugation PDB is unavailable.
    """
    source = working_dir / SOLVATED_CONJUGATE_PDB
    alias = working_dir / SOLVATED_SYSTEM_PDB
    if not source.exists():
        raise FileNotFoundError(
            f"Conjugation PDB '{source}' is required to create '{alias.name}'. "
            "Run 'polyzymd build' without --skip-build or restore the conjugation PDB."
        )

    if not alias.exists() or not filecmp.cmp(source, alias, shallow=False):
        alias.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, alias)
    return alias


def resolve_skip_build_artifacts(sim_config: object, working_dir: Path) -> tuple[Path, Path]:
    """Validate and resolve pre-built OpenMM simulation inputs.

    Parameters
    ----------
    sim_config : object
        Simulation configuration used to determine route-specific validation.
    working_dir : Path
        Replicate working directory to inspect.

    Returns
    -------
    tuple[Path, Path]
        ``(pdb_path, system_xml_path)`` for OpenMM startup.

    Raises
    ------
    FileNotFoundError
        If required handoff files are missing.
    """
    system_path = working_dir / "system.xml"
    if not system_path.exists():
        raise FileNotFoundError(
            f"Pre-built OpenMM system '{system_path}' not found. "
            "Run 'polyzymd build' first or remove --skip-build."
        )

    if conjugation_enabled(sim_config):
        return ensure_conjugation_pdb_alias(working_dir), system_path

    pdb_path = working_dir / SOLVATED_SYSTEM_PDB
    if pdb_path.exists():
        return pdb_path, system_path

    raise FileNotFoundError(
        f"Pre-built topology '{pdb_path}' not found. Run 'polyzymd build' first or remove "
        "--skip-build."
    )


def write_openmm_system_xml(
    *,
    builder: Any,
    sim_config: object,
    working_dir: Path,
    apply_configured_restraints: bool = True,
) -> Path:
    """Serialize an OpenMM ``System`` from a prepared builder.

    Parameters
    ----------
    builder : Any
        Builder exposing ``get_openmm_components()``.
    sim_config : object
        Simulation configuration containing optional restraints.
    working_dir : Path
        Directory where ``system.xml`` should be written.

    Returns
    -------
    Path
        Path to the serialized ``system.xml`` file.
    """
    omm_topology, omm_system, _omm_positions = builder.get_openmm_components()

    restraints = getattr(sim_config, "restraints", None) if apply_configured_restraints else None
    if restraints:
        from polyzymd.core.restraints import RestraintFactory, apply_restraints

        restraint_defs = []
        for restraint in restraints:
            if not restraint.enabled:
                continue
            restraint_defs.append(RestraintFactory.from_config(restraint.model_dump()))
        if restraint_defs:
            apply_restraints(restraint_defs, omm_topology, omm_system)

    from openmm import XmlSerializer

    working_dir.mkdir(parents=True, exist_ok=True)
    system_path = working_dir / "system.xml"
    with open(system_path, "w") as f:
        f.write(XmlSerializer.serialize(omm_system))
    return system_path


def _write_solute_audit(
    *,
    builder: Any,
    solute_dir: Path,
    pdb_path: Path,
    system_path: Path,
    exact_export_bundle: Any | None = None,
) -> Path:
    """Write deterministic evidence for an isolated-primary OpenMM build."""
    from openmm import XmlSerializer, unit

    system = XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    force_classes = Counter(type(system.getForce(i)).__name__ for i in range(system.getNumForces()))
    restraint_force_count = 0
    for index in range(system.getNumForces()):
        force = system.getForce(index)
        force_class = type(force).__name__
        tagged = str(force.getName()).startswith(_RESTRAINT_FORCE_PREFIX)
        legacy_custom_bond = force_class == "CustomBondForce" and (
            force.getEnergyFunction() in _RESTRAINT_CUSTOM_BOND_EXPRESSIONS
        )
        restraint_force_count += int(
            tagged or force_class == "CustomExternalForce" or legacy_custom_bond
        )
    if restraint_force_count:
        raise RuntimeError(
            "Solute-scope parameterization unexpectedly contains "
            f"{restraint_force_count} restraint-like force(s)"
        )
    periodic = any(
        bool(system.getForce(i).usesPeriodicBoundaryConditions())
        for i in range(system.getNumForces())
    )
    box_vectors_nm = None
    if periodic:
        vectors = system.getDefaultPeriodicBoxVectors()
        box_vectors_nm = [
            [float(vector[index].value_in_unit(unit.nanometer)) for index in range(3)]
            for vector in vectors
        ]
    atom_lines = [
        line
        for line in pdb_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    pdb_atom_count = len(atom_lines)
    particle_count = int(system.getNumParticles())
    if pdb_atom_count != particle_count:
        raise RuntimeError(
            "Solute coordinate atom count does not match serialized OpenMM System particle "
            f"count: PDB={pdb_atom_count}, System={particle_count}"
        )
    coordinates = [
        tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
        for line in atom_lines
    ]
    if not all(math.isfinite(value) for coordinate in coordinates for value in coordinate):
        raise RuntimeError("Solute coordinate artifact contains non-finite coordinates")
    if exact_export_bundle is not None:
        topology_atoms = int(exact_export_bundle.topology.getNumAtoms())
        sidecar = exact_export_bundle.sidecar
        if topology_atoms != particle_count or sidecar.particle_count != particle_count:
            raise RuntimeError("Exact solute topology, sidecar, and System counts do not match")
        pdb_identities = tuple(
            (line[12:16].strip(), line[17:20].strip(), line[21].strip(), line[22:26].strip())
            for line in atom_lines
        )
        exact_identities = tuple(
            (atom.name, atom.residue_name, atom.chain_id, atom.residue_id) for atom in sidecar.atoms
        )
        if pdb_identities != exact_identities:
            raise RuntimeError("Exact solute PDB atom identities/order do not match the sidecar")
        # Re-validation proves the stored order/topology/parameter hashes are exact.
        type(sidecar).model_validate(sidecar.model_dump())

    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    residues = {(line[17:20].strip().upper(), line[21:27]) for line in atom_lines}
    component_counts = {
        "substrate": int(getattr(builder, "_n_substrate_molecules", 0)),
        "free_polymers": int(getattr(builder, "_n_polymer_chains", 0)),
        "water": sum(name in WATER_RESIDUES for name, _identity in residues),
        "ions": sum(name in COMMON_ION_RESIDUES for name, _identity in residues),
        "liquids": sum(name in COMMON_SOLVENT_RESIDUES for name, _identity in residues),
    }
    if any(component_counts.values()):
        raise RuntimeError(
            "Solute scope contains excluded secondary components: " f"{component_counts}"
        )
    audit = {
        "scope": "solute",
        "route": "native_openmm_glycam" if exact_export_bundle is not None else "generic_openff",
        "supported": True,
        "pdb_path": pdb_path.name,
        "system_path": system_path.name,
        "pdb_atom_count": pdb_atom_count,
        "system_particle_count": particle_count,
        "component_counts": component_counts,
        "force_classes": dict(sorted(force_classes.items())),
        "barostat_count": sum(count for name, count in force_classes.items() if "Barostat" in name),
        "restraint_force_count": restraint_force_count,
        "periodic": periodic,
        "box_vectors_nm": box_vectors_nm,
        "sha256": {"solute.pdb": digest(pdb_path), "system.xml": digest(system_path)},
    }
    if exact_export_bundle is not None:
        if not periodic:
            raise RuntimeError("Exact native solute must retain its periodic PME box")
        audit["exact_hashes"] = {
            "atom_order": sidecar.atom_order_hash,
            "topology": sidecar.topology_hash,
            "particles": sidecar.particle_hash,
            "exceptions": sidecar.exception_hash,
        }
        audit["preparation_only_warning"] = exact_export_bundle.audit.get(
            "preparation_only_warning"
        )
    if audit["barostat_count"] != 0:
        raise RuntimeError("Solute-scope parameterization unexpectedly contains a barostat")
    audit_path = solute_dir / "openmm_build_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return audit_path


def _reject_disconnected_hetatm(pdb_path: Path) -> None:
    """Fail when any HETATM atom is outside the primary covalent graph."""
    lines = pdb_path.read_text(encoding="utf-8").splitlines()
    atom_serials = {
        int(line[6:11]) for line in lines if line.startswith("ATOM  ") and line[6:11].strip()
    }
    hetatm_serials = {
        int(line[6:11]) for line in lines if line.startswith("HETATM") and line[6:11].strip()
    }
    if not hetatm_serials:
        return
    adjacency = parse_pdb_conect_adjacency(pdb_path)
    reachable = set(atom_serials)
    pending = list(atom_serials)
    while pending:
        for neighbor in adjacency.get(pending.pop(), ()):
            if neighbor not in reachable:
                reachable.add(neighbor)
                pending.append(neighbor)
    disconnected = sorted(hetatm_serials - reachable)
    if disconnected:
        raise ValueError(
            "build scope 'solute' found disconnected/unowned HETATM atoms in the assembled "
            f"primary component (serials: {disconnected}); supply a primary-structure-only input"
        )


def build_openmm_artifacts(
    *,
    sim_config: object,
    working_dir: Path,
    polymer_seed: int,
    write_system: bool = True,
    scope: BuildScope = BuildScope.SYSTEM,
) -> OpenMMBuildArtifacts:
    """Build OpenMM handoff artifacts through the standard or conjugation route.

    Parameters
    ----------
    sim_config : object
        Validated simulation configuration.
    working_dir : Path
        Replicate working directory for build outputs.
    polymer_seed : int
        Replicate-specific seed used for polymer placement.
    write_system : bool, optional
        If ``True``, serialize ``system.xml``. Conjugation builds refresh the
        standard PDB alias whenever the route-specific PDB exists.

    Returns
    -------
    OpenMMBuildArtifacts
        Build metadata and standardized OpenMM handoff paths.
    """
    if scope is BuildScope.STRUCTURE and not conjugation_enabled(sim_config):
        raise ValueError("build scope 'structure' requires conjugation.enabled: true")

    if conjugation_enabled(sim_config):
        from polyzymd.builders.conjugation import build_conjugate_from_config
        from polyzymd.builders.conjugation.native_openmm_glycam import native_glycam_enabled
        from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings

        workflow_settings = ConjugatedPolymerSystemSettings(create_final_interchange=True)
        if scope is not BuildScope.SYSTEM:
            workflow_settings = workflow_settings.model_copy(update={"build_scope": scope})
        if native_glycam_enabled(sim_config):
            workflow_settings = workflow_settings.model_copy(
                update={"pdb_fragment_output_mode": "experimental_pablo"}
            )
        if scope is BuildScope.SOLUTE:
            with tempfile.TemporaryDirectory(prefix="polyzymd-solute-") as temporary_dir:
                result = build_conjugate_from_config(
                    sim_config,
                    output_dir=Path(temporary_dir),
                    settings=workflow_settings,
                    free_polymer_seed=polymer_seed,
                )
        else:
            result = build_conjugate_from_config(
                sim_config,
                output_dir=working_dir,
                settings=workflow_settings,
                free_polymer_seed=polymer_seed,
            )
        builder = result.system_builder
        if scope is BuildScope.STRUCTURE:
            pdb_path = result.crosslinked_conjugate_pdb_path
            if pdb_path is None or not Path(pdb_path).is_file():
                raise RuntimeError(
                    "Structure-scope conjugation did not produce assembled_crosslinked.pdb"
                )
            return OpenMMBuildArtifacts(
                builder=None,
                interchange=None,
                result=result,
                pdb_path=Path(pdb_path),
                system_path=None,
                conjugation_enabled=True,
                scope=scope,
            )
        if scope is BuildScope.SOLUTE:
            exact_bundle = getattr(result, "exact_export_bundle", None)
            if builder is None:
                raise RuntimeError("Solute-scope route did not return a SystemBuilder")
            solute_dir = working_dir / "solute"
            solute_dir.mkdir(parents=True, exist_ok=True)
            pdb_path = solute_dir / "solute.pdb"
            if exact_bundle is not None:
                exact_bundle.write_pdb(pdb_path)
            else:
                builder.save_topology(pdb_path)
            _reject_disconnected_hetatm(pdb_path)
            if exact_bundle is not None:
                from openmm import XmlSerializer

                system_path = solute_dir / "system.xml"
                system_path.write_text(
                    XmlSerializer.serialize(exact_bundle.to_openmm()), encoding="utf-8"
                )
                exact_bundle.sidecar_path = exact_bundle.sidecar.save(
                    solute_dir / "exact_openmm_exceptions.json"
                )
                exact_bundle.audit_path = solute_dir / "native_openmm_glycam_audit.json"
                exact_bundle.audit_path.write_text(
                    json.dumps(exact_bundle.audit, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8",
                )
            else:
                system_path = write_openmm_system_xml(
                    builder=builder,
                    sim_config=sim_config,
                    working_dir=solute_dir,
                    apply_configured_restraints=False,
                )
            _write_solute_audit(
                builder=builder,
                solute_dir=solute_dir,
                pdb_path=pdb_path,
                system_path=system_path,
                exact_export_bundle=exact_bundle,
            )
            return OpenMMBuildArtifacts(
                builder=builder,
                interchange=None if exact_bundle is not None else builder.interchange,
                result=result,
                exact_export_bundle=exact_bundle,
                pdb_path=pdb_path,
                system_path=system_path,
                conjugation_enabled=True,
                scope=scope,
            )
        if builder is not None:
            if write_system:
                system_path = write_openmm_system_xml(
                    builder=builder,
                    sim_config=sim_config,
                    working_dir=working_dir,
                )
            else:
                system_path = working_dir / "system.xml"
        else:
            result.require_final_interchange()
            system_path = working_dir / "system.xml"

        conjugation_pdb_path = working_dir / SOLVATED_CONJUGATE_PDB
        pdb_path = (
            ensure_conjugation_pdb_alias(working_dir)
            if write_system or conjugation_pdb_path.exists()
            else working_dir / SOLVATED_SYSTEM_PDB
        )
        return OpenMMBuildArtifacts(
            builder=builder,
            interchange=None,
            result=result,
            exact_export_bundle=getattr(result, "exact_export_bundle", None),
            pdb_path=pdb_path,
            system_path=system_path,
            conjugation_enabled=True,
            scope=scope,
        )

    from polyzymd.builders.system_builder import SystemBuilder

    builder = SystemBuilder.from_config(sim_config)
    if scope is BuildScope.SOLUTE:
        builder.build_isolated_primary_from_config(sim_config)
        solute_dir = working_dir / "solute"
        solute_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = solute_dir / "solute.pdb"
        builder.save_topology(pdb_path)
        _reject_disconnected_hetatm(pdb_path)
        system_path = write_openmm_system_xml(
            builder=builder,
            sim_config=sim_config,
            working_dir=solute_dir,
            apply_configured_restraints=False,
        )
        _write_solute_audit(
            builder=builder, solute_dir=solute_dir, pdb_path=pdb_path, system_path=system_path
        )
        return OpenMMBuildArtifacts(
            builder=builder,
            interchange=builder.interchange,
            result=None,
            exact_export_bundle=None,
            pdb_path=pdb_path,
            system_path=system_path,
            conjugation_enabled=False,
            scope=scope,
        )
    interchange = builder.build_from_config(
        config=sim_config,
        working_dir=working_dir,
        polymer_seed=polymer_seed,
    )
    if write_system:
        system_path = write_openmm_system_xml(
            builder=builder,
            sim_config=sim_config,
            working_dir=working_dir,
        )
    else:
        system_path = working_dir / "system.xml"
    return OpenMMBuildArtifacts(
        builder=builder,
        interchange=interchange,
        result=None,
        exact_export_bundle=None,
        pdb_path=working_dir / SOLVATED_SYSTEM_PDB,
        system_path=system_path,
        conjugation_enabled=False,
        scope=scope,
    )
