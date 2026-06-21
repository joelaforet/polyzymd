"""Helpers for building and validating OpenMM handoff artifacts."""

from __future__ import annotations

import filecmp
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SOLVATED_SYSTEM_PDB = "solvated_system.pdb"
SOLVATED_CONJUGATE_PDB = "solvated_conjugate_free_polymers.pdb"


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
    system_path: Path
    conjugation_enabled: bool

    def get_component_info(self) -> Any:
        """Return component metadata from the route-specific build object.

        Returns
        -------
        Any
            Component metadata suitable for exporter calls.
        """
        if self.result is not None:
            return self.result.get_component_info()
        return self.builder.get_component_info()

    def require_final_interchange(self) -> Any:
        """Return the final Interchange from the route-specific build object.

        Returns
        -------
        Any
            Final Interchange object for export.
        """
        if self.result is not None:
            return self.result.require_final_interchange()
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


def write_openmm_system_xml(*, builder: Any, sim_config: object, working_dir: Path) -> Path:
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

    restraints = getattr(sim_config, "restraints", None)
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


def build_openmm_artifacts(
    *,
    sim_config: object,
    working_dir: Path,
    polymer_seed: int,
    write_system: bool = True,
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
    if conjugation_enabled(sim_config):
        from polyzymd.builders.conjugation import build_conjugate_from_config
        from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings

        result = build_conjugate_from_config(
            sim_config,
            output_dir=working_dir,
            settings=ConjugatedPolymerSystemSettings(create_final_interchange=True),
            free_polymer_seed=polymer_seed,
        )
        builder = result.system_builder
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
            pdb_path=pdb_path,
            system_path=system_path,
            conjugation_enabled=True,
        )

    from polyzymd.builders.system_builder import SystemBuilder

    builder = SystemBuilder.from_config(sim_config)
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
        pdb_path=working_dir / SOLVATED_SYSTEM_PDB,
        system_path=system_path,
        conjugation_enabled=False,
    )
