"""Opt-in integrated physics smoke test for conjugation construction.

This is the acceptance test for the real construction path: Polymerist recipe
generation, Polymerist PDB parsing, Packmol placement, crosslinked PDB writing,
Pablo ingestion with explicit LYX/NHX crosslink configuration, OpenFF
Interchange parameterization, and restrained vacuum OpenMM minimization plus MD.

Run on Blanca GPU resources with a command similar to::

    module load slurm/blanca
    salloc ...
    PYTHONNOUSERSITE=1 POLYZYMD_RUN_CONJUGATION_PABLO_SMOKE=1 \
      pixi run -e sim-cuda-12-4 pytest \
      tests/test_conjugation_integrated_smoke.py -v

When the opt-in environment variable is set, scientific or software blockers
after workflow start are reported as hard failures with the artifact directory.
"""

from __future__ import annotations

import os
import shutil
import time
import traceback
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import TypeVar

import pytest

from polyzymd.builders.conjugation._assembly import (
    PackmolModifierPlacementSettings,
    place_modifier_with_packmol,
)
from polyzymd.builders.conjugation._linkage import (
    NhsLysModifierLinker,
    PabloCrosslinkRequirement,
    require_pablo_crosslink_requirement,
)
from polyzymd.builders.conjugation._relaxation import (
    VacuumSmokeSettings,
    run_restrained_vacuum_smoke,
)
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestionResult, PabloIngestor
from polyzymd.builders.conjugation.pablo.parameterization import (
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.polymer.polymerist import generated_fragment_from_polymerist_pdb
from polyzymd.builders.conjugation.polymer.recipe import generate_multi_residue_molecule
from polyzymd.builders.conjugation.structure.pdb import (
    CrosslinkedPdbAssemblyOptions,
    PdbAtomRecord,
    write_crosslinked_pdb,
)
from tests._support.conjugation_polymer_recipes import sbma_egpma_nhs_recipe

T = TypeVar("T")

SMOKE_ENV_VAR = "POLYZYMD_RUN_CONJUGATION_PABLO_SMOKE"
SMOKE_PLATFORM_ENV_VAR = "POLYZYMD_CONJUGATION_SMOKE_PLATFORM"
SMOKE_MINIMIZATION_ITERS_ENV_VAR = "POLYZYMD_CONJUGATION_SMOKE_MIN_ITERS"
SMOKE_NVT_STEPS_ENV_VAR = "POLYZYMD_CONJUGATION_SMOKE_NVT_STEPS"
SMOKE_POLYMERIST_RETRIES_ENV_VAR = "POLYZYMD_CONJUGATION_SMOKE_POLYMERIST_RETRIES"
POC_PROTEIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "polyzymd"
    / "builders"
    / "conjugation"
    / "poc"
    / "data"
    / "NH3_terminal_His_proton_updated.pdb"
)
TARGET_LYSINE_RESIDUE = 23
ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


@pytest.mark.slow
@pytest.mark.conjugation_stack
def test_opt_in_integrated_conjugation_physics_smoke(tmp_path: Path):
    """Run the real integrated conjugation construction and vacuum MD smoke."""
    platform_name = _require_conjugation_stack_or_skip()
    artifact_dir = tmp_path / "conjugation-integrated-smoke"
    artifact_dir.mkdir(parents=True, exist_ok=True)

    if not POC_PROTEIN_PATH.exists():
        pytest.fail(f"Required lightweight POC protein fixture is missing: {POC_PROTEIN_PATH}")

    recipe = sbma_egpma_nhs_recipe(length=3, seed=51, reactive_monomer_index=1)
    generation = _run_stage(
        "SBMA/EGPMA/NHS multi-residue generation",
        artifact_dir,
        lambda: generate_multi_residue_molecule(
            recipe,
            artifact_dir / "polymerist-cache",
            force_regenerate=True,
            max_retries=_int_env(SMOKE_POLYMERIST_RETRIES_ENV_VAR, default=1),
        ),
    )
    if generation.pdb_path is None:
        _fail_blocker(
            "SBMA/EGPMA/NHS multi-residue generation",
            RuntimeError("Generation backend did not return a generated PDB path"),
            artifact_dir,
        )

    modifier = _run_stage(
        "Polymerist PDB parsing",
        artifact_dir,
        lambda: generated_fragment_from_polymerist_pdb(
            generation.pdb_path,
            recipe=recipe,
            sequence=generation.sequence,
            name="smoke-sbma-egpma-nhs",
        ),
    )
    linker = NhsLysModifierLinker(target_residue_number=TARGET_LYSINE_RESIDUE)
    resolved_plan = linker.resolve_plan(POC_PROTEIN_PATH, modifier)
    pablo_policy = _explicit_lyx_nhx_policy(resolved_plan.pablo_crosslink_requirement)

    crosslink_validation = _run_stage(
        "Explicit LYX/NHX Pablo crosslink validation",
        artifact_dir,
        lambda: require_pablo_crosslink_requirement(
            pablo_policy,
            resolved_plan.pablo_crosslink_requirement,
        ),
    )
    assert crosslink_validation.residues == ("LYX", "NHX")

    placement = _run_stage(
        "Packmol modifier placement",
        artifact_dir,
        lambda: place_modifier_with_packmol(
            POC_PROTEIN_PATH,
            modifier,
            linker,
            artifact_dir,
            settings=PackmolModifierPlacementSettings(nloop=500),
        ),
    )
    assert placement.packmol_output_path.exists()

    crosslinked_pdb = artifact_dir / "assembled_crosslinked.pdb"
    assembly = _run_stage(
        "Crosslinked PDB assembly",
        artifact_dir,
        lambda: write_crosslinked_pdb(
            POC_PROTEIN_PATH,
            placement.placed_modifier,
            linker.attachment(POC_PROTEIN_PATH),
            crosslinked_pdb,
            CrosslinkedPdbAssemblyOptions(),
        ),
    )
    assert assembly.output_path == crosslinked_pdb
    assert crosslinked_pdb.exists()

    pablo_result = _run_stage(
        "Pablo LYX/NHX crosslinked PDB ingestion",
        artifact_dir,
        lambda: PabloIngestor(policy=pablo_policy).ingest_structure(
            crosslinked_pdb,
            output_dir=artifact_dir,
        ),
    )
    if not pablo_result.success or pablo_result.topology is None:
        _fail_blocker(
            "Pablo LYX/NHX crosslinked PDB ingestion",
            RuntimeError(_pablo_failure_summary(pablo_result)),
            artifact_dir,
        )

    parameterization = _run_stage(
        "OpenFF Interchange parameterization",
        artifact_dir,
        lambda: create_interchange_from_pablo_topology(pablo_result.topology),
    )
    if not parameterization.success or parameterization.interchange is None:
        _fail_blocker(
            "OpenFF Interchange parameterization",
            RuntimeError("OpenFF Interchange did not return a parameterized system"),
            artifact_dir,
        )

    restrained_indices = _protein_heavy_atom_indices(crosslinked_pdb)
    smoke_result = _run_stage(
        "Restrained vacuum OpenMM minimization and MD",
        artifact_dir,
        lambda: run_restrained_vacuum_smoke(
            parameterization.interchange,
            artifact_dir,
            protein_heavy_atom_indices=restrained_indices,
            settings=VacuumSmokeSettings(
                minimization_max_iterations=_int_env(
                    SMOKE_MINIMIZATION_ITERS_ENV_VAR,
                    default=5,
                ),
                nvt_steps=_int_env(SMOKE_NVT_STEPS_ENV_VAR, default=2),
                platform_name=platform_name,
            ),
        ),
    )

    assert smoke_result.success
    assert smoke_result.nvt_steps == _int_env(SMOKE_NVT_STEPS_ENV_VAR, default=2)
    assert smoke_result.smoke_json_path.exists()
    assert smoke_result.minimized_pdb_path is not None and smoke_result.minimized_pdb_path.exists()
    assert (
        smoke_result.equilibrated_pdb_path is not None
        and smoke_result.equilibrated_pdb_path.exists()
    )


def _require_conjugation_stack_or_skip() -> str:
    """Skip unless the opt-in chemistry stack and OpenMM platform are available."""
    if os.environ.get(SMOKE_ENV_VAR) != "1":
        pytest.skip(f"Set {SMOKE_ENV_VAR}=1 to run the Pablo conjugation physics smoke")
    if shutil.which("packmol") is None:
        pytest.skip("Packmol binary is not available on PATH")

    pytest.importorskip("polymerist", exc_type=ImportError)
    pytest.importorskip("openff.pablo")
    pytest.importorskip("openff.toolkit")
    pytest.importorskip("openff.interchange")
    pytest.importorskip("openmm")

    return _select_usable_openmm_platform_or_skip()


def _select_usable_openmm_platform_or_skip() -> str:
    """Return a usable OpenMM platform, preferring CUDA unless overridden."""
    import openmm
    from openmm import unit

    available_platforms = {
        openmm.Platform.getPlatform(index).getName()
        for index in range(openmm.Platform.getNumPlatforms())
    }
    requested = os.environ.get(SMOKE_PLATFORM_ENV_VAR)
    candidates = (requested,) if requested else ("CUDA", "OpenCL", "CPU", "Reference")
    errors: list[str] = []
    for platform_name in candidates:
        if platform_name is None or platform_name not in available_platforms:
            continue
        try:
            _validate_openmm_platform_context(openmm, unit, platform_name)
        except Exception as exc:  # noqa: BLE001 - OpenMM platform failures vary by resource
            errors.append(f"{platform_name}: {exc}")
            continue
        return platform_name

    details = "; ".join(errors) if errors else "no requested platform is registered"
    pytest.skip(f"No usable OpenMM platform for conjugation smoke: {details}")


def _validate_openmm_platform_context(openmm: object, unit: object, platform_name: str) -> None:
    """Create a one-particle context to prove the platform resource is usable."""
    system = openmm.System()
    system.addParticle(39.9)
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName(platform_name)
    context = openmm.Context(system, integrator, platform)
    context.setPositions([[0.0, 0.0, 0.0]] * unit.nanometer)
    context.getState(getEnergy=True)
    del context
    del integrator


def _explicit_lyx_nhx_policy(requirement: PabloCrosslinkRequirement) -> SimpleNamespace:
    """Build an explicit Pablo policy for the realized LYX/NHX linkage."""
    return SimpleNamespace(
        lookup_policy="auto_download",
        ccd_cache_directory=None,
        residue_definition_files=(),
        use_canonical_atom_names=False,
        crosslinks=[
            SimpleNamespace(
                residues=requirement.residues,
                linking_atoms=requirement.linking_atoms,
                leaving_atoms=requirement.leaving_atoms,
                bond_order=requirement.bond_order,
            )
        ],
    )


def _run_stage(stage: str, artifact_dir: Path, func: Callable[[], T]) -> T:
    """Run one real workflow stage and fail with an actionable artifact pointer."""
    start = time.monotonic()
    print(f"[conjugation-smoke] START {stage}", flush=True)
    try:
        result = func()
    except Exception as exc:  # noqa: BLE001 - third-party stack errors need normalized context
        _fail_blocker(stage, exc, artifact_dir)
    elapsed = time.monotonic() - start
    print(f"[conjugation-smoke] DONE  {stage} ({elapsed:.1f} s)", flush=True)
    return result


def _int_env(name: str, *, default: int) -> int:
    """Read a positive integer environment setting for the opt-in smoke."""
    value = os.environ.get(name)
    if value is None or value == "":
        return default
    parsed = int(value)
    if parsed < 0:
        raise ValueError(f"{name} must be non-negative, got {parsed}")
    return parsed


def _fail_blocker(stage: str, exc: Exception, artifact_dir: Path) -> None:
    """Fail the opt-in smoke with explicit blocker guidance and saved artifacts."""
    artifact_dir.mkdir(parents=True, exist_ok=True)
    report_path = artifact_dir / "conjugation_smoke_failure.txt"
    report_path.write_text(
        f"Stage: {stage}\n"
        f"Artifact directory: {artifact_dir}\n"
        f"Exception type: {type(exc).__name__}\n"
        f"Exception: {exc}\n\n"
        f"Traceback:\n{traceback.format_exc()}\n",
        encoding="utf-8",
    )
    pytest.fail(
        f"Integrated conjugation physics smoke failed during {stage}. "
        f"{_blocker_guidance(stage, exc)} Artifacts were left under {artifact_dir}. "
        f"Failure report: {report_path}. Original error: {type(exc).__name__}: {exc}"
    )


def _blocker_guidance(stage: str, exc: Exception) -> str:
    """Return stage-specific guidance for expected current blockers."""
    message = str(exc).lower()
    if "pablo" in stage.lower() or any(token in message for token in ("lyx", "nhx", "ccd")):
        return (
            "This means Pablo could not ingest the LYX/NHX crosslinked PDB. Add or fix "
            "Pablo residue definitions for LYX and NHX, then verify the explicit "
            "ccd_pablo.crosslinks entry names both product residues and linking atoms."
        )
    if "interchange" in stage.lower() or any(
        token in message for token in ("force field", "parameter", "template")
    ):
        return (
            "This means the Pablo topology reached parameterization, but OpenFF force field "
            "coverage is incomplete. Add templates, charges, and bonded parameters for the "
            "LYX/NHX product residues before claiming a full physics smoke."
        )
    if "openmm" in stage.lower() or "vacuum" in stage.lower():
        return (
            "This means a parameterized system was built, but restrained vacuum EM/MD did "
            "not complete with finite states. Inspect the minimized and smoke JSON artifacts."
        )
    return "This is a real workflow blocker, not a dependency guard skip."


def _pablo_failure_summary(result: PabloIngestionResult) -> str:
    """Return a compact Pablo failure summary from diagnostics."""
    diagnostics = [f"{diag.code}: {diag.message}" for diag in result.diagnostics]
    joined = "; ".join(diagnostics) if diagnostics else "no diagnostics were reported"
    return f"Pablo returned success={result.success}, topology={result.topology}: {joined}"


def _protein_heavy_atom_indices(pdb_path: Path) -> tuple[int, ...]:
    """Return zero-based chain-A non-hydrogen atom indices from a PDB artifact."""
    indices: list[int] = []
    atom_index = 0
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(ATOM_RECORD_PREFIXES):
                continue
            atom = PdbAtomRecord.from_pdb_line(line, atom_index=atom_index)
            if atom.chain_id == "A" and not _is_hydrogen(atom):
                indices.append(atom_index)
            atom_index += 1
    if not indices:
        pytest.fail(f"No chain-A heavy atoms were found for restraints in {pdb_path}")
    return tuple(indices)


def _is_hydrogen(atom: PdbAtomRecord) -> bool:
    """Return whether a parsed atom record is hydrogen."""
    return (atom.element or "").upper() == "H" or atom.atom_name.upper().startswith("H")
