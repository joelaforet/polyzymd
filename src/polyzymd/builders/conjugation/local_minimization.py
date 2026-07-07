"""Post-crosslink local minimization helpers for conjugate products."""

from __future__ import annotations

import json
import math
import traceback
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import PabloCrosslinkRequirement, parse_pdb_atom_records
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor
from polyzymd.builders.conjugation.pablo.parameterization import (
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.relaxation._diagnostics import (
    _positions_to_numpy,
    validate_finite_energy,
    validate_finite_positions,
)
from polyzymd.builders.conjugation.relaxation._openmm_system import (
    _add_positional_restraints,
    _openmm_positions_from_interchange,
    _state_energy_kj_mol,
)
from polyzymd.builders.conjugation.structure.parsing import parse_pdb_conect_records
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.config.schema import ConjugationCcdCrosslinkConfig, ConjugationCcdPabloPolicyConfig


class LocalLinkageAtomSelector(BaseModel):
    """Fixed PDB atom selector for one product-state linkage atom."""

    serial: int | None = Field(None, ge=1)
    chain_id: str
    residue_name: str
    residue_number: int
    atom_name: str


class LocalLinkageSelectors(BaseModel):
    """Selectors and target geometry for one product-state conjugation linkage."""

    protein_link_selector: LocalLinkageAtomSelector
    modifier_link_selector: LocalLinkageAtomSelector
    target_bond_length_angstrom: float = Field(1.33, gt=0)
    label: str = "linkage"


class LocalLinkageGeometryMetrics(BaseModel):
    """Geometry and connectivity metrics for one product-state linkage."""

    label: str
    protein_modifier_distance_angstrom: float
    target_bond_length_angstrom: float
    distance_error_angstrom: float
    reciprocal_product_link_conect: bool
    passes: bool
    failures: tuple[str, ...] = Field(default_factory=tuple)


class LocalGeometryMetrics(BaseModel):
    """Geometry and connectivity metrics for all product-state linkages."""

    linkages: tuple[LocalLinkageGeometryMetrics, ...]
    passes: bool
    failures: tuple[str, ...] = Field(default_factory=tuple)


class LocalMinimizationSettings(BaseModel):
    """Settings for restrained post-crosslink OpenMM/OpenFF/Pablo minimization."""

    linkages: tuple[LocalLinkageSelectors, ...] = Field(default_factory=tuple)
    polymer_mobile_residue_window: int = Field(1, ge=0)
    restraint_k_kj_mol_nm2: float = Field(1_000_000.0, gt=0)
    minimization_tolerance_kj_mol_nm: float = Field(1.0, gt=0)
    minimization_max_iterations: int = Field(500, ge=0)
    max_protein_displacement_angstrom: float = Field(0.05, ge=0)
    max_linkage_distance_error_angstrom: float = Field(0.35, ge=0)
    platform_name: str | None = None
    relaxed_pdb_name: str = "crosslinked_relaxed.pdb"
    result_json_name: str = "local_minimization_result.json"
    use_default_pablo_crosslink: bool = True


class LocalMinimizationResult(BaseModel):
    """Result record for a local minimization attempt."""

    success: bool
    input_pdb: Path
    output_dir: Path
    result_json_path: Path | None = None
    relaxed_pdb_path: Path | None = None
    before_geometry: LocalGeometryMetrics
    after_geometry: LocalGeometryMetrics | None = None
    platform_name: str | None = None
    mobile_atom_count: int = 0
    restrained_atom_count: int = 0
    energy_before_min_kj_mol: float | None = None
    energy_after_min_kj_mol: float | None = None
    max_restrained_protein_displacement_angstrom: float | None = None
    mean_restrained_protein_displacement_angstrom: float | None = None
    blocker: str | None = None
    blocker_traceback: str | None = None
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    def save(self, path: Path | str | None = None) -> Path:
        """Write this result as JSON and return the output path."""
        output_path = (
            Path(path) if path is not None else self.output_dir / "local_minimization_result.json"
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.result_json_path = output_path
        output_path.write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")
        return output_path


def analyze_crosslink_geometry(
    pdb_path: Path | str,
    *,
    settings: LocalMinimizationSettings | None = None,
) -> LocalGeometryMetrics:
    """Measure product-link distances and required reciprocal ``CONECT`` records."""
    minimization_settings = settings or LocalMinimizationSettings()
    if not minimization_settings.linkages:
        raise ValueError("Local minimization geometry analysis requires at least one linkage")
    atoms = parse_pdb_atom_records(pdb_path)
    conect = _parse_conect_records(Path(pdb_path))
    linkage_metrics: list[LocalLinkageGeometryMetrics] = []
    all_failures: list[str] = []
    for linkage in minimization_settings.linkages:
        protein_atom = _select_atom(atoms, linkage.protein_link_selector)
        modifier_atom = _select_atom(atoms, linkage.modifier_link_selector)
        distance = _distance(_atom_position(protein_atom), _atom_position(modifier_atom))
        error = abs(distance - linkage.target_bond_length_angstrom)
        reciprocal = _has_conect(conect, protein_atom.serial, modifier_atom.serial) and _has_conect(
            conect, modifier_atom.serial, protein_atom.serial
        )
        failures: list[str] = []
        if error > minimization_settings.max_linkage_distance_error_angstrom:
            failures.append(
                f"{linkage.label} distance error {error:.3f} A exceeds "
                f"{minimization_settings.max_linkage_distance_error_angstrom:.3f} A"
            )
        if not reciprocal:
            failures.append(f"{linkage.label} product link CONECT is not reciprocal")
        all_failures.extend(failures)
        linkage_metrics.append(
            LocalLinkageGeometryMetrics(
                label=linkage.label,
                protein_modifier_distance_angstrom=distance,
                target_bond_length_angstrom=linkage.target_bond_length_angstrom,
                distance_error_angstrom=error,
                reciprocal_product_link_conect=reciprocal,
                passes=not failures,
                failures=tuple(failures),
            )
        )

    return LocalGeometryMetrics(
        linkages=tuple(linkage_metrics),
        passes=not all_failures,
        failures=tuple(all_failures),
    )


def run_post_crosslink_local_minimization(
    pdb_path: Path | str,
    output_dir: Path | str | None = None,
    *,
    settings: LocalMinimizationSettings | None = None,
    pablo_policy: Any | None = None,
    pablo_crosslink_requirement: Any | None = None,
    product_state_pablo_library: Any | None = None,
    resolved_plan: Any | None = None,
) -> LocalMinimizationResult:
    """Run restrained local minimization through OpenMM/OpenFF/Pablo.

    The function never falls back to RDKit/UFF. If Pablo ingestion or OpenFF
    parameterization fails, the returned result records the blocker and traceback.
    """
    minimization_settings = settings or LocalMinimizationSettings()
    input_pdb = Path(pdb_path)
    artifact_dir = Path(output_dir) if output_dir is not None else input_pdb.parent
    artifact_dir.mkdir(parents=True, exist_ok=True)
    before_geometry = analyze_crosslink_geometry(input_pdb, settings=minimization_settings)

    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as openmm_unit
    except ImportError as exc:
        return _blocked_result(
            input_pdb,
            artifact_dir,
            before_geometry,
            exc,
            "OpenMM is required for local minimization; run in the build pixi environment.",
            minimization_settings,
        )

    try:
        custom_residue_library = _extract_residue_library(product_state_pablo_library)
        policy = pablo_policy
        if policy is None and custom_residue_library is None:
            policy = _product_state_pablo_policy(
                input_pdb,
                minimization_settings,
                pablo_crosslink_requirement=pablo_crosslink_requirement,
                resolved_plan=resolved_plan,
            )
        ingestion = PabloIngestor(policy).ingest_structure(
            input_pdb,
            output_dir=artifact_dir,
            residue_library=custom_residue_library,
        )
        if not ingestion.success or ingestion.topology is None:
            raise RuntimeError(_pablo_failure_message(ingestion))

        interchange_result = create_interchange_from_pablo_topology(
            ingestion.topology,
            charge_from_molecules=(),
        )
        interchange = interchange_result.interchange
        topology = interchange.to_openmm_topology()
        initial_positions = _openmm_positions_from_interchange(interchange, openmm_unit)
        initial_positions_array = _positions_to_numpy(initial_positions, openmm_unit)

        atoms = parse_pdb_atom_records(input_pdb)
        if len(atoms) != len(initial_positions_array):
            raise RuntimeError(
                "Pablo/OpenFF atom count does not match the input PDB order: "
                f"PDB has {len(atoms)} atoms but positions have {len(initial_positions_array)}"
            )

        mobile_indices = _mobile_atom_indices(atoms, minimization_settings)
        restrained_indices = _restrained_protein_indices(atoms, mobile_indices)
        if not restrained_indices:
            raise RuntimeError("No nonparticipating protein atoms were selected for restraints")

        system = interchange.to_openmm_system()
        _add_positional_restraints(
            system,
            initial_positions,
            restrained_indices,
            minimization_settings.restraint_k_kj_mol_nm2,
            openmm,
            openmm_unit,
        )
        simulation, platform_name = _build_simulation_with_platform_fallback(
            openmm,
            openmm_app,
            topology,
            system,
            minimization_settings,
            openmm_unit,
        )
        simulation.context.setPositions(initial_positions)

        energy_before = _state_energy_kj_mol(
            simulation.context.getState(getEnergy=True), openmm_unit
        )
        validate_finite_energy(energy_before, label="energy_before_min_kj_mol")
        openmm.LocalEnergyMinimizer.minimize(
            simulation.context,
            tolerance=minimization_settings.minimization_tolerance_kj_mol_nm
            * openmm_unit.kilojoule_per_mole
            / openmm_unit.nanometer,
            maxIterations=minimization_settings.minimization_max_iterations,
        )
        state_after = simulation.context.getState(getEnergy=True, getPositions=True)
        energy_after = _state_energy_kj_mol(state_after, openmm_unit)
        minimized_positions = state_after.getPositions(asNumpy=True)
        validate_finite_energy(energy_after, label="energy_after_min_kj_mol")
        validate_finite_positions(minimized_positions, openmm_unit, label="minimized_positions")

        minimized_angstrom = _positions_to_numpy(minimized_positions, openmm_unit) * 10.0
        relaxed_path = artifact_dir / minimization_settings.relaxed_pdb_name
        write_pdb_with_replaced_coordinates(input_pdb, minimized_angstrom, relaxed_path)
        after_geometry = analyze_crosslink_geometry(relaxed_path, settings=minimization_settings)
        max_displacement, mean_displacement = _protein_restraint_displacements_angstrom(
            initial_positions_array * 10.0,
            minimized_angstrom,
            restrained_indices,
        )

        result = LocalMinimizationResult(
            success=after_geometry.passes
            and (max_displacement <= minimization_settings.max_protein_displacement_angstrom),
            input_pdb=input_pdb,
            output_dir=artifact_dir,
            relaxed_pdb_path=relaxed_path,
            before_geometry=before_geometry,
            after_geometry=after_geometry,
            platform_name=platform_name,
            mobile_atom_count=len(mobile_indices),
            restrained_atom_count=len(restrained_indices),
            energy_before_min_kj_mol=energy_before,
            energy_after_min_kj_mol=energy_after,
            max_restrained_protein_displacement_angstrom=max_displacement,
            mean_restrained_protein_displacement_angstrom=mean_displacement,
            diagnostics=tuple(_success_diagnostics(after_geometry, max_displacement)),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        result = _blocked_result(
            input_pdb,
            artifact_dir,
            before_geometry,
            exc,
            "OpenMM/OpenFF/Pablo local minimization could not be completed.",
            minimization_settings,
        )

    result.save(artifact_dir / minimization_settings.result_json_name)
    return result


def _build_simulation_with_platform_fallback(
    openmm: Any,
    openmm_app: Any,
    topology: Any,
    system: Any,
    settings: LocalMinimizationSettings,
    openmm_unit: Any,
) -> tuple[Any, str]:
    """Create an OpenMM Simulation, retrying CPU/Reference after accelerator failures."""
    names = (
        (settings.platform_name,)
        if settings.platform_name is not None
        else ("CUDA", "OpenCL", "CPU", "Reference")
    )
    errors: list[str] = []
    for name in names:
        if name is None:
            continue
        try:
            platform = openmm.Platform.getPlatformByName(name)
            integrator = openmm.VerletIntegrator(0.001 * openmm_unit.picoseconds)
            simulation = openmm_app.Simulation(topology, system, integrator, platform)
            platform_name = platform.getName() if hasattr(platform, "getName") else str(platform)
            return simulation, platform_name
        except Exception as exc:  # noqa: BLE001 - OpenMM platform/context errors vary
            errors.append(f"{name}: {exc}")
            if settings.platform_name is not None:
                break
    raise RuntimeError(
        "No suitable OpenMM platform found for local minimization: " + "; ".join(errors)
    )


def write_pdb_with_replaced_coordinates(
    template_pdb_path: Path | str,
    coordinates_angstrom: Any,
    output_path: Path | str,
) -> Path:
    """Write a PDB by replacing only ATOM/HETATM coordinate columns."""
    template_path = Path(template_pdb_path)
    output = Path(output_path)
    coords = np.asarray(coordinates_angstrom, dtype=float)
    atom_count = sum(
        1
        for line in template_path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )
    if coords.shape != (atom_count, 3):
        raise ValueError(
            "Coordinate replacement requires one XYZ triplet per PDB atom: "
            f"expected {(atom_count, 3)}, got {coords.shape}"
        )
    if not np.all(np.isfinite(coords)):
        raise ValueError("Replacement coordinates contain non-finite values")

    output.parent.mkdir(parents=True, exist_ok=True)
    coord_index = 0
    out_lines: list[str] = []
    with template_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(("ATOM", "HETATM")):
                x_coord, y_coord, z_coord = coords[coord_index]
                padded = f"{line:<80}"
                line = f"{padded[:30]}{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}" f"{padded[54:]}"
                coord_index += 1
            out_lines.append(line)
    output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return output


def product_state_pablo_crosslink_requirement(
    product_pdb_path: Path | str,
    *,
    settings: LocalMinimizationSettings | None = None,
    pablo_crosslink_requirement: Any | None = None,
    resolved_plan: Any | None = None,
    leaving_atoms: tuple[tuple[str, ...], tuple[str, ...]] = ((), ()),
) -> PabloCrosslinkRequirement:
    """Return a Pablo crosslink requirement for an already-modified product PDB.

    PolyzyMD removes reactant-side leaving atoms before writing the local
    minimization input PDB. Pablo should therefore be told to match the emitted
    product residues and crosslink atom names, not to remove reactant-state atoms
    such as lysine NZ hydrogens again.
    """
    minimization_settings = settings or LocalMinimizationSettings()
    source = _coerce_pablo_crosslink_requirement(
        pablo_crosslink_requirement or resolved_plan,
        default_settings=minimization_settings,
    )
    requirement = PabloCrosslinkRequirement(
        residues=source.residues,
        linking_atoms=source.linking_atoms,
        leaving_atoms=leaving_atoms,
        bond_order=source.bond_order,
    )
    _validate_product_state_crosslink_atoms(product_pdb_path, requirement)
    return requirement


def build_product_state_pablo_policy(
    product_pdb_path: Path | str,
    *,
    settings: LocalMinimizationSettings | None = None,
    pablo_crosslink_requirement: Any | None = None,
    resolved_plan: Any | None = None,
) -> ConjugationCcdPabloPolicyConfig | None:
    """Build a Pablo policy for product-state local minimization ingestion."""
    minimization_settings = settings or LocalMinimizationSettings()
    if not minimization_settings.use_default_pablo_crosslink and not (
        pablo_crosslink_requirement or resolved_plan
    ):
        return None
    requirement = product_state_pablo_crosslink_requirement(
        product_pdb_path,
        settings=minimization_settings,
        pablo_crosslink_requirement=pablo_crosslink_requirement,
        resolved_plan=resolved_plan,
    )
    return _policy_from_crosslink_requirement(requirement)


def _blocked_result(
    input_pdb: Path,
    output_dir: Path,
    before_geometry: LocalGeometryMetrics,
    exc: BaseException,
    message: str,
    settings: LocalMinimizationSettings,
) -> LocalMinimizationResult:
    """Build a failed result preserving the full exception traceback."""
    return LocalMinimizationResult(
        success=False,
        input_pdb=input_pdb,
        output_dir=output_dir,
        before_geometry=before_geometry,
        blocker=f"{message} {type(exc).__name__}: {exc}",
        blocker_traceback="".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
        diagnostics=(
            "No RDKit/UFF fallback was attempted; minimization requires the OpenMM/OpenFF/Pablo route.",
            "The input PDB metadata and CONECT records were preserved; only geometry was measured.",
            "Product-state Pablo crosslinks use empty leaving atom groups only when no custom product library is supplied.",
            f"Pablo crosslink policy enabled: {settings.use_default_pablo_crosslink}",
        ),
    )


def _extract_residue_library(product_state_pablo_library: Any | None) -> Any | None:
    """Return a Pablo cache from a product-state library wrapper or raw cache."""
    if product_state_pablo_library is None:
        return None
    return getattr(product_state_pablo_library, "residue_library", product_state_pablo_library)


def _product_state_pablo_policy(
    product_pdb_path: Path | str,
    settings: LocalMinimizationSettings,
    *,
    pablo_crosslink_requirement: Any | None = None,
    resolved_plan: Any | None = None,
) -> ConjugationCcdPabloPolicyConfig | None:
    """Return the Pablo policy used by the POC local minimizer."""
    return build_product_state_pablo_policy(
        product_pdb_path,
        settings=settings,
        pablo_crosslink_requirement=pablo_crosslink_requirement,
        resolved_plan=resolved_plan,
    )


def _policy_from_crosslink_requirement(
    requirement: PabloCrosslinkRequirement,
) -> ConjugationCcdPabloPolicyConfig:
    """Build a Pablo policy from a normalized crosslink requirement."""
    return ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            ConjugationCcdCrosslinkConfig(
                residues=requirement.residues,
                linking_atoms=requirement.linking_atoms,
                leaving_atoms=requirement.leaving_atoms,
                bond_order=requirement.bond_order,
            )
        ]
    )


def _coerce_pablo_crosslink_requirement(
    value: Any | None,
    *,
    default_settings: LocalMinimizationSettings,
) -> PabloCrosslinkRequirement:
    """Normalize resolved plans, requirements, dicts, or objects to a requirement."""
    if value is None:
        if not default_settings.linkages:
            raise ValueError("Cannot infer a Pablo crosslink without resolved linkage settings")
        linkage = default_settings.linkages[0]
        return PabloCrosslinkRequirement(
            residues=(
                linkage.protein_link_selector.residue_name,
                linkage.modifier_link_selector.residue_name,
            ),
            linking_atoms=(
                linkage.protein_link_selector.atom_name,
                linkage.modifier_link_selector.atom_name,
            ),
            leaving_atoms=((), ()),
            bond_order=1,
        )
    if isinstance(value, PabloCrosslinkRequirement):
        return value
    plan_requirement = getattr(value, "pablo_crosslink_requirement", None)
    if plan_requirement is not None:
        return _coerce_pablo_crosslink_requirement(
            plan_requirement,
            default_settings=default_settings,
        )
    if isinstance(value, dict):
        return PabloCrosslinkRequirement.model_validate(value)
    return PabloCrosslinkRequirement(
        residues=value.residues,
        linking_atoms=value.linking_atoms,
        leaving_atoms=getattr(value, "leaving_atoms", ((), ())),
        bond_order=getattr(value, "bond_order", 1),
    )


def _validate_product_state_crosslink_atoms(
    product_pdb_path: Path | str,
    requirement: PabloCrosslinkRequirement,
) -> None:
    """Validate that product-state residue and crosslink atom names exist in the PDB."""
    atoms = parse_pdb_atom_records(product_pdb_path)
    missing: list[str] = []
    for residue_name, atom_name in zip(requirement.residues, requirement.linking_atoms):
        if not any(
            atom.residue_name.upper() == residue_name and atom.atom_name.upper() == atom_name
            for atom in atoms
        ):
            missing.append(f"{residue_name}:{atom_name}")
    if missing:
        raise ValueError(
            "Product-state Pablo crosslink atoms were not found in the emitted PDB: "
            + ", ".join(missing)
        )


def _pablo_failure_message(ingestion: Any) -> str:
    """Return a compact Pablo failure message from an ingestion result."""
    errors = [
        diagnostic
        for diagnostic in getattr(ingestion, "diagnostics", ())
        if getattr(getattr(diagnostic, "severity", None), "value", None) == "error"
    ]
    if not errors:
        return "Pablo ingestion failed without a detailed diagnostic"
    details = getattr(errors[-1], "details", {})
    return str(details.get("error") or getattr(errors[-1], "message", "Pablo ingestion failed"))


def _success_diagnostics(
    geometry: LocalGeometryMetrics,
    max_displacement_angstrom: float,
) -> list[str]:
    """Return human-readable success/failure diagnostics after minimization."""
    diagnostics = [
        "OpenMM/OpenFF/Pablo local minimization completed without RDKit/UFF fallback",
        f"Maximum restrained protein displacement was {max_displacement_angstrom:.4f} A",
    ]
    if geometry.passes:
        diagnostics.append("Relaxed geometry passed product-link validation")
    else:
        diagnostics.extend(geometry.failures)
    return diagnostics


def _select_atom(
    atoms: tuple[PdbAtomRecord, ...], selector: LocalLinkageAtomSelector
) -> PdbAtomRecord:
    """Return exactly one atom matching a fixed selector."""
    matches = [atom for atom in atoms if _selector_matches(atom, selector)]
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one atom for selector "
            f"{selector.model_dump(mode='json')}, found {len(matches)}"
        )
    return matches[0]


def _selector_matches(atom: PdbAtomRecord, selector: LocalLinkageAtomSelector) -> bool:
    """Return whether an atom matches a crosslink selector."""
    return (
        (selector.serial is None or atom.serial == selector.serial)
        and atom.chain_id.upper() == selector.chain_id.upper()
        and atom.residue_name.upper() == selector.residue_name.upper()
        and atom.residue_number == selector.residue_number
        and atom.atom_name.upper() == selector.atom_name.upper()
    )


def _atom_position(atom: PdbAtomRecord) -> np.ndarray:
    """Return an atom position in Angstrom."""
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def _distance(first: np.ndarray, second: np.ndarray) -> float:
    """Return a Cartesian distance in Angstrom."""
    return float(np.linalg.norm(first - second))


def _angle_degrees(first: np.ndarray, vertex: np.ndarray, third: np.ndarray) -> float:
    """Return the angle first-vertex-third in degrees."""
    vector_a = first - vertex
    vector_b = third - vertex
    norm_product = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
    if norm_product == 0.0:
        raise ValueError("Cannot compute angle with coincident atoms")
    cosine = float(np.dot(vector_a, vector_b) / norm_product)
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def _parse_conect_records(path: Path) -> dict[int, set[int]]:
    """Parse PDB ``CONECT`` records into a serial adjacency map."""
    return {source: set(targets) for source, targets in parse_pdb_conect_records(path).items()}


def _has_conect(conect: dict[int, set[int]], source: int | None, target: int | None) -> bool:
    """Return whether a directed ``CONECT`` entry exists."""
    if source is None or target is None:
        return False
    return target in conect.get(source, set())


def _mobile_atom_indices(
    atoms: tuple[PdbAtomRecord, ...],
    settings: LocalMinimizationSettings,
) -> tuple[int, ...]:
    """Select protein link residues and local modifier residues allowed to move."""
    mobile: set[int] = set()
    for linkage in settings.linkages:
        protein_selector = linkage.protein_link_selector
        modifier_selector = linkage.modifier_link_selector
        modifier_residue_min = (
            modifier_selector.residue_number - settings.polymer_mobile_residue_window
        )
        modifier_residue_max = (
            modifier_selector.residue_number + settings.polymer_mobile_residue_window
        )
        for atom in atoms:
            if atom.atom_index is None:
                continue
            if _same_selector_residue(atom, protein_selector):
                mobile.add(atom.atom_index)
            if (
                atom.chain_id.upper() == modifier_selector.chain_id.upper()
                and modifier_residue_min <= atom.residue_number <= modifier_residue_max
            ):
                mobile.add(atom.atom_index)
    return tuple(sorted(mobile))


def _same_selector_residue(atom: PdbAtomRecord, selector: LocalLinkageAtomSelector) -> bool:
    """Return whether an atom is in the selected residue."""
    return (
        atom.chain_id.upper() == selector.chain_id.upper()
        and atom.residue_name.upper() == selector.residue_name.upper()
        and atom.residue_number == selector.residue_number
    )


def _restrained_protein_indices(
    atoms: tuple[PdbAtomRecord, ...],
    mobile_indices: tuple[int, ...],
) -> tuple[int, ...]:
    """Select all non-mobile protein atoms for strong positional restraints."""
    mobile = set(mobile_indices)
    indices = [
        atom.atom_index
        for atom in atoms
        if atom.atom_index is not None
        and atom.chain_id.upper() == "A"
        and atom.atom_index not in mobile
    ]
    return tuple(sorted(indices))


def _protein_restraint_displacements_angstrom(
    initial_angstrom: np.ndarray,
    final_angstrom: np.ndarray,
    restrained_indices: tuple[int, ...],
) -> tuple[float, float]:
    """Return max and mean displacement for restrained protein atoms."""
    if not restrained_indices:
        return 0.0, 0.0
    displacement = np.linalg.norm(
        final_angstrom[list(restrained_indices)] - initial_angstrom[list(restrained_indices)],
        axis=1,
    )
    return float(np.max(displacement)), float(np.mean(displacement))
