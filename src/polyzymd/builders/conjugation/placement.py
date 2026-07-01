"""Packmol placement helpers for conjugation workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import (
    ModifierLinker,
    ProteinLinkSite,
    ResolvedAttachmentPlan,
)
from polyzymd.builders.conjugation.polymer import GeneratedPolymerFragment
from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_records as parse_structure_pdb_atom_records,
)
from polyzymd.builders.conjugation.structure.parsing import pdb_coordinates
from polyzymd.builders.conjugation.structure.pdb import (
    PdbAtomRecord,
    PlacedPolymerFragment,
)

_SHIFT_PADDING_ANGSTROM = 10.0
_BOX_PADDING_ANGSTROM = 30.0


class PackmolModifierPlacementSettings(BaseModel):
    """Settings for constrained modifier placement with Packmol."""

    tolerance_angstrom: float = Field(2.0, gt=0)
    reactive_sphere_radius_angstrom: float = Field(5.0, gt=0)
    target_bond_length_angstrom: float = Field(1.33, gt=0)
    movebadrandom: bool = True
    nloop: int = Field(500, gt=0)
    work_dir_name: str = "packmol_modifier_placement"


class PackmolModifierPlacementResult(BaseModel):
    """Placed modifier plus Packmol diagnostics and artifacts."""

    placed_modifier: PlacedPolymerFragment
    packmol_input_path: Path
    packmol_output_path: Path
    protein_sterics_pdb_path: Path
    modifier_pdb_path: Path
    packmol_input_text: str
    target_bond_length_angstrom: float
    placed_bond_length_angstrom: float
    min_modifier_protein_distance_angstrom: float
    excluded_protein_atom_names: tuple[str, ...] = Field(default_factory=tuple)


def place_modifier_with_packmol(
    protein_pdb_path: Path | str,
    modifier: GeneratedPolymerFragment,
    linker: ModifierLinker,
    output_dir: Path | str,
    *,
    settings: PackmolModifierPlacementSettings | None = None,
    run_packmol_func: Any | None = None,
) -> PackmolModifierPlacementResult:
    """Place a generated modifier near a protein reactive atom using Packmol.

    Packmol receives the protein as a fixed steric obstacle and the modifier as
    a movable structure. The modifier reactive atom is constrained inside a
    sphere centered on the protein linker atom, then the whole modifier is
    translated so the reactive atoms are exactly at the requested bond length.

    Parameters
    ----------
    protein_pdb_path : pathlib.Path or str
        Protein PDB containing the target linker atom.
    modifier : GeneratedPolymerFragment
        Generated modifier or polymer fragment before placement.
    linker : ModifierLinker
        Linker strategy that resolves the protein site.
    output_dir : pathlib.Path or str
        Directory for Packmol artifacts.
    settings : PackmolModifierPlacementSettings or None, optional
        Placement settings, by default ``None``.
    run_packmol_func : callable or None, optional
        Optional Packmol executor for tests, by default ``None``.

    Returns
    -------
    PackmolModifierPlacementResult
        Placed modifier and artifact paths.
    """
    plan_builder = getattr(linker, "resolve_plan", None)
    if not callable(plan_builder):
        raise TypeError("Modifier linker must provide resolve_plan() for Packmol placement")
    plan = plan_builder(protein_pdb_path, modifier)
    if settings is not None:
        plan = plan.model_copy(
            update={"target_bond_length_angstrom": settings.target_bond_length_angstrom}
        )
    return place_modifier_with_resolved_plan(
        protein_pdb_path,
        modifier,
        plan,
        output_dir,
        settings=settings,
        run_packmol_func=run_packmol_func,
    )


def place_modifier_with_resolved_plan(
    protein_pdb_path: Path | str,
    modifier: GeneratedPolymerFragment,
    plan: ResolvedAttachmentPlan,
    output_dir: Path | str,
    *,
    settings: PackmolModifierPlacementSettings | None = None,
    run_packmol_func: Any | None = None,
) -> PackmolModifierPlacementResult:
    """Place a generated modifier using a resolved explicit linkage plan.

    Parameters
    ----------
    protein_pdb_path : pathlib.Path or str
        Protein PDB containing the resolved protein link atom.
    modifier : GeneratedPolymerFragment
        Generated modifier or polymer fragment before placement.
    plan : ResolvedAttachmentPlan
        Resolved atom-level protein/modifier linkage plan.
    output_dir : pathlib.Path or str
        Directory for Packmol artifacts.
    settings : PackmolModifierPlacementSettings or None, optional
        Placement settings, by default ``None``.
    run_packmol_func : callable or None, optional
        Optional Packmol executor for tests, by default ``None``.

    Returns
    -------
    PackmolModifierPlacementResult
        Placed modifier and artifact paths.
    """
    from polyzymd.utils.packmol import build_packmol_input, run_packmol

    placement_settings = settings or PackmolModifierPlacementSettings()
    protein_path = Path(protein_pdb_path).resolve()
    work_dir = (Path(output_dir) / placement_settings.work_dir_name).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    protein_atoms = _parse_pdb_atoms(protein_path)
    retained_modifier_atoms = _retained_modifier_atoms(modifier, plan=plan)
    reactive_atom = plan.modifier_link_atom

    protein_steric_atoms, excluded_names = _protein_steric_atoms_from_plan(protein_atoms, plan)
    protein_coords = _coords_from_atoms(protein_steric_atoms)
    modifier_retained_coords = _coords_from_atoms(retained_modifier_atoms)
    modifier_full_coords = _coords_from_atoms(tuple(atom.to_pdb_atom() for atom in modifier.atoms))

    coord_shift = _coordinate_shift(protein_coords, modifier_retained_coords)
    shifted_protein = protein_coords + coord_shift
    shifted_modifier = modifier_retained_coords + coord_shift
    shifted_site = _coord(plan.protein_link_atom) + coord_shift
    box_size = np.vstack((shifted_protein, shifted_modifier)).max(axis=0) + _BOX_PADDING_ANGSTROM

    protein_sterics_path = work_dir / "protein_fixed_sterics.pdb"
    modifier_path = work_dir / "modifier_retained.pdb"
    _write_simple_pdb(protein_sterics_path, shifted_protein, ["C"] * len(shifted_protein))
    _write_simple_pdb(
        modifier_path,
        shifted_modifier,
        [(atom.element or "C").strip() or "C" for atom in retained_modifier_atoms],
    )

    reactive_local_index = _retained_local_index(retained_modifier_atoms, reactive_atom)
    structure_extra_lines = [
        [
            f"atoms {reactive_local_index + 1}",
            (
                "inside sphere "
                f"{shifted_site[0]:.6f} {shifted_site[1]:.6f} {shifted_site[2]:.6f} "
                f"{placement_settings.reactive_sphere_radius_angstrom:.1f}"
            ),
            "end atoms",
        ]
    ]

    packmol_input_text = build_packmol_input(
        molecule_pdb_paths=[str(modifier_path)],
        molecule_counts=[1],
        box_size_angstrom=box_size,
        tolerance_angstrom=placement_settings.tolerance_angstrom,
        solute_pdb_path=str(protein_sterics_path),
        use_pbc=False,
        movebadrandom=placement_settings.movebadrandom,
        nloop=placement_settings.nloop,
        structure_extra_lines=structure_extra_lines,
    )
    packmol_input_path = work_dir / "packmol.inp"
    packmol_input_path.write_text(packmol_input_text, encoding="utf-8")

    executor = run_packmol_func or run_packmol
    packmol_output_path = Path(executor(packmol_input_text, work_dir))
    packed_coords = _read_pdb_coords(packmol_output_path)
    expected_atoms = len(shifted_protein) + len(retained_modifier_atoms)
    if len(packed_coords) != expected_atoms:
        raise RuntimeError(
            f"Packmol output has {len(packed_coords)} atoms, expected {expected_atoms} "
            f"(protein={len(shifted_protein)}, modifier={len(retained_modifier_atoms)})"
        )

    packed_modifier_coords = packed_coords[len(shifted_protein) :] - coord_shift
    transformed_full_coords = _transform_full_modifier(
        modifier_retained_coords,
        packed_modifier_coords,
        modifier_full_coords,
    )
    transformed_full_coords = _snap_reactive_atom_to_bond_length(
        transformed_full_coords,
        reactive_atom,
        plan.protein_link_atom,
        plan.target_bond_length_angstrom,
    )
    placed_fragment = _placed_fragment_from_coords(modifier, transformed_full_coords)

    placed_reactive = _coord(resolve_modifier_reactive_atom_from_placed(placed_fragment))
    placed_bond_length = float(np.linalg.norm(placed_reactive - _coord(plan.protein_link_atom)))
    min_distance = _minimum_distance(
        _coords_from_atoms(tuple(placed_fragment.atoms)),
        _coords_from_atoms(protein_steric_atoms),
    )

    return PackmolModifierPlacementResult(
        placed_modifier=placed_fragment,
        packmol_input_path=packmol_input_path,
        packmol_output_path=packmol_output_path,
        protein_sterics_pdb_path=protein_sterics_path,
        modifier_pdb_path=modifier_path,
        packmol_input_text=packmol_input_text,
        target_bond_length_angstrom=plan.target_bond_length_angstrom,
        placed_bond_length_angstrom=placed_bond_length,
        min_modifier_protein_distance_angstrom=min_distance,
        excluded_protein_atom_names=excluded_names,
    )


def place_modifiers_with_resolved_plans(
    protein_pdb_path: Path | str,
    modifiers: tuple[GeneratedPolymerFragment, ...],
    plans: tuple[ResolvedAttachmentPlan, ...],
    output_dir: Path | str,
    *,
    settings: PackmolModifierPlacementSettings | None = None,
    run_packmol_func: Any | None = None,
) -> tuple[PackmolModifierPlacementResult, ...]:
    """Place multiple generated modifiers with one joint Packmol run.

    Packmol receives one fixed protein sterics structure and one movable
    retained-fragment structure per attachment. All fragment restraints are
    solved together so protein-fragment and fragment-fragment sterics are
    optimized in the same packing problem.
    """
    from polyzymd.utils.packmol import build_packmol_input, run_packmol

    if not modifiers or len(modifiers) != len(plans):
        raise ValueError("Joint Packmol placement requires aligned modifiers and plans")

    placement_settings = settings or PackmolModifierPlacementSettings()
    protein_path = Path(protein_pdb_path).resolve()
    artifact_dir = Path(output_dir).resolve()
    work_dir = (artifact_dir / placement_settings.work_dir_name).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    protein_atoms = _parse_pdb_atoms(protein_path)
    protein_steric_atoms, excluded_names = _protein_steric_atoms_from_plans(
        protein_atoms,
        plans,
    )
    protein_coords = _coords_from_atoms(protein_steric_atoms)

    retained_atom_groups = tuple(
        _retained_modifier_atoms(modifier, plan=plan)
        for modifier, plan in zip(modifiers, plans, strict=True)
    )
    retained_coord_groups = tuple(_coords_from_atoms(atoms) for atoms in retained_atom_groups)
    full_coord_groups = tuple(
        _coords_from_atoms(tuple(atom.to_pdb_atom() for atom in modifier.atoms))
        for modifier in modifiers
    )

    coord_shift = _coordinate_shift(protein_coords, *retained_coord_groups)
    shifted_protein = protein_coords + coord_shift
    shifted_retained_groups = tuple(coords + coord_shift for coords in retained_coord_groups)
    box_size = (
        np.vstack((shifted_protein, *shifted_retained_groups)).max(axis=0) + _BOX_PADDING_ANGSTROM
    )

    protein_sterics_path = work_dir / "protein_fixed_sterics.pdb"
    _write_simple_pdb(protein_sterics_path, shifted_protein, ["C"] * len(shifted_protein))

    modifier_paths: list[Path] = []
    structure_extra_lines: list[list[str]] = []
    for index, (retained_atoms, shifted_modifier, plan) in enumerate(
        zip(retained_atom_groups, shifted_retained_groups, plans, strict=True),
        start=1,
    ):
        fragment_work_dir = (
            artifact_dir / f"placement_{index:02d}" / placement_settings.work_dir_name
        ).resolve()
        fragment_work_dir.mkdir(parents=True, exist_ok=True)
        modifier_path = fragment_work_dir / "modifier_retained.pdb"
        _write_simple_pdb(
            modifier_path,
            shifted_modifier,
            [(atom.element or "C").strip() or "C" for atom in retained_atoms],
        )
        modifier_paths.append(modifier_path)

        reactive_local_index = _retained_local_index(retained_atoms, plan.modifier_link_atom)
        shifted_site = _coord(plan.protein_link_atom) + coord_shift
        structure_extra_lines.append(
            [
                f"atoms {reactive_local_index + 1}",
                (
                    "inside sphere "
                    f"{shifted_site[0]:.6f} {shifted_site[1]:.6f} {shifted_site[2]:.6f} "
                    f"{placement_settings.reactive_sphere_radius_angstrom:.1f}"
                ),
                "end atoms",
            ]
        )

    packmol_input_text = build_packmol_input(
        molecule_pdb_paths=[str(path) for path in modifier_paths],
        molecule_counts=[1] * len(modifier_paths),
        box_size_angstrom=box_size,
        tolerance_angstrom=placement_settings.tolerance_angstrom,
        solute_pdb_path=str(protein_sterics_path),
        use_pbc=False,
        movebadrandom=placement_settings.movebadrandom,
        nloop=placement_settings.nloop,
        structure_extra_lines=structure_extra_lines,
    )
    packmol_input_path = work_dir / "packmol.inp"
    packmol_input_path.write_text(packmol_input_text, encoding="utf-8")

    executor = run_packmol_func or run_packmol
    packmol_output_path = Path(executor(packmol_input_text, work_dir))
    packed_coords = _read_pdb_coords(packmol_output_path)
    retained_counts = tuple(len(atoms) for atoms in retained_atom_groups)
    expected_atoms = len(shifted_protein) + sum(retained_counts)
    if len(packed_coords) != expected_atoms:
        raise RuntimeError(
            f"Packmol output has {len(packed_coords)} atoms, expected {expected_atoms} "
            f"(protein={len(shifted_protein)}, modifiers={retained_counts})"
        )

    results: list[PackmolModifierPlacementResult] = []
    start = len(shifted_protein)
    for modifier, plan, retained_coords, full_coords, count, modifier_path in zip(
        modifiers,
        plans,
        retained_coord_groups,
        full_coord_groups,
        retained_counts,
        modifier_paths,
        strict=True,
    ):
        packed_modifier_coords = packed_coords[start : start + count] - coord_shift
        start += count

        transformed_full_coords = _transform_full_modifier(
            retained_coords,
            packed_modifier_coords,
            full_coords,
        )
        transformed_full_coords = _snap_reactive_atom_to_bond_length(
            transformed_full_coords,
            plan.modifier_link_atom,
            plan.protein_link_atom,
            plan.target_bond_length_angstrom,
        )
        placed_fragment = _placed_fragment_from_coords(modifier, transformed_full_coords)

        placed_reactive = _coord(resolve_modifier_reactive_atom_from_placed(placed_fragment))
        placed_bond_length = float(np.linalg.norm(placed_reactive - _coord(plan.protein_link_atom)))
        min_distance = _minimum_distance(
            _coords_from_atoms(tuple(placed_fragment.atoms)),
            _coords_from_atoms(protein_steric_atoms),
        )
        results.append(
            PackmolModifierPlacementResult(
                placed_modifier=placed_fragment,
                packmol_input_path=packmol_input_path,
                packmol_output_path=packmol_output_path,
                protein_sterics_pdb_path=protein_sterics_path,
                modifier_pdb_path=modifier_path,
                packmol_input_text=packmol_input_text,
                target_bond_length_angstrom=plan.target_bond_length_angstrom,
                placed_bond_length_angstrom=placed_bond_length,
                min_modifier_protein_distance_angstrom=min_distance,
                excluded_protein_atom_names=excluded_names,
            )
        )

    return tuple(results)


def resolve_modifier_reactive_atom_from_placed(fragment: PlacedPolymerFragment) -> PdbAtomRecord:
    """Resolve a reactive atom from a placed fragment."""
    matches = [
        atom
        for atom in fragment.atoms
        if (
            fragment.reactive_atom_serial is not None
            and atom.serial == fragment.reactive_atom_serial
        )
        or (
            fragment.reactive_atom_index is not None
            and atom.atom_index == fragment.reactive_atom_index
        )
        or (
            fragment.reactive_atom_name is not None
            and atom.atom_name.upper() == fragment.reactive_atom_name.upper()
        )
    ]
    if len(matches) != 1:
        raise ValueError(f"Placed modifier reactive atom selector resolved {len(matches)} atoms")
    return matches[0]


def _parse_pdb_atoms(path: Path) -> tuple[PdbAtomRecord, ...]:
    """Parse atom records from a PDB file."""
    return parse_structure_pdb_atom_records(path, require_atoms=True)


def _protein_steric_atoms(
    atoms: tuple[PdbAtomRecord, ...], site: ProteinLinkSite
) -> tuple[tuple[PdbAtomRecord, ...], tuple[str, ...]]:
    """Return fixed protein sterics after excluding bond-forming atoms."""
    excluded_indices = {site.atom.atom_index}
    excluded_indices.update(atom.atom_index for atom in site.removable_hydrogens)
    kept = tuple(atom for atom in atoms if atom.atom_index not in excluded_indices)
    excluded_names = tuple(
        atom.atom_name
        for atom in atoms
        if atom.atom_index in excluded_indices and atom.atom_index is not None
    )
    return kept, excluded_names


def _protein_steric_atoms_from_plan(
    atoms: tuple[PdbAtomRecord, ...], plan: ResolvedAttachmentPlan
) -> tuple[tuple[PdbAtomRecord, ...], tuple[str, ...]]:
    """Return fixed protein sterics after excluding resolved linkage atoms."""
    excluded_identities = {_atom_identity(plan.protein_link_atom)}
    excluded_identities.update(_atom_identity(atom) for atom in plan.protein_leaving_atoms)
    kept = tuple(atom for atom in atoms if _atom_identity(atom) not in excluded_identities)
    excluded_names = tuple(
        atom.atom_name for atom in atoms if _atom_identity(atom) in excluded_identities
    )
    return kept, excluded_names


def _protein_steric_atoms_from_plans(
    atoms: tuple[PdbAtomRecord, ...], plans: tuple[ResolvedAttachmentPlan, ...]
) -> tuple[tuple[PdbAtomRecord, ...], tuple[str, ...]]:
    """Return fixed protein sterics after excluding all resolved linkage atoms."""
    excluded_identities = set()
    for plan in plans:
        excluded_identities.add(_atom_identity(plan.protein_link_atom))
        excluded_identities.update(_atom_identity(atom) for atom in plan.protein_leaving_atoms)
    kept = tuple(atom for atom in atoms if _atom_identity(atom) not in excluded_identities)
    excluded_names = tuple(
        atom.atom_name for atom in atoms if _atom_identity(atom) in excluded_identities
    )
    return kept, excluded_names


def _retained_modifier_atoms(
    modifier: GeneratedPolymerFragment,
    *,
    plan: ResolvedAttachmentPlan | None = None,
) -> tuple[PdbAtomRecord, ...]:
    """Return modifier atoms retained for Packmol placement."""
    placed = modifier.to_placed_fragment()
    if plan is not None:
        leaving_identities = {_atom_identity(atom) for atom in plan.modifier_leaving_atoms}
        return tuple(
            atom for atom in placed.atoms if _atom_identity(atom) not in leaving_identities
        )

    leaving_serials = set(placed.leaving_atom_serials)
    leaving_indices = set(placed.leaving_atom_indices)
    leaving_names = {name.upper() for name in placed.leaving_atom_names}
    return tuple(
        atom
        for atom in placed.atoms
        if not (
            (atom.serial is not None and atom.serial in leaving_serials)
            or atom.atom_index in leaving_indices
            or atom.atom_name.upper() in leaving_names
        )
    )


def _coords_from_atoms(atoms: tuple[PdbAtomRecord, ...]) -> np.ndarray:
    """Return atom coordinates as an ``(N, 3)`` array in Angstrom."""
    return np.asarray([_coord(atom) for atom in atoms], dtype=float)


def _coord(atom: PdbAtomRecord) -> np.ndarray:
    """Return one atom coordinate vector."""
    return np.asarray([atom.x, atom.y, atom.z], dtype=float)


def _coordinate_shift(*coord_arrays: np.ndarray) -> np.ndarray:
    """Return a positive-octant coordinate shift for Packmol."""
    combined = np.vstack(coord_arrays)
    return -combined.min(axis=0) + _SHIFT_PADDING_ANGSTROM


def _write_simple_pdb(path: Path, coords: np.ndarray, elements: list[str]) -> None:
    """Write a minimal PDB file for Packmol steric placement."""
    lines = []
    for index, (xyz, element) in enumerate(zip(coords, elements, strict=True), start=1):
        elem = (element or "C").strip().upper()[:2] or "C"
        name = f"{elem}{index:>3d}"[:4]
        lines.append(
            f"HETATM{index:5d} {name:<4s} UNK A   1    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {elem:>2s}\n"
        )
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def _read_pdb_coords(path: Path) -> np.ndarray:
    """Read coordinates from PDB atom records."""
    return np.asarray(pdb_coordinates(path, require_atoms=False), dtype=float)


def _retained_local_index(
    retained_atoms: tuple[PdbAtomRecord, ...], reactive_atom: PdbAtomRecord
) -> int:
    """Return the Packmol-local index of the reactive atom."""
    for index, atom in enumerate(retained_atoms):
        if _same_atom(atom, reactive_atom):
            return index
    raise ValueError("Reactive atom is not among retained modifier atoms")


def _transform_full_modifier(
    retained_source: np.ndarray,
    retained_target: np.ndarray,
    full_source: np.ndarray,
) -> np.ndarray:
    """Apply the retained-atom rigid transform to all modifier atoms."""
    if len(retained_source) < 3:
        translation = retained_target.mean(axis=0) - retained_source.mean(axis=0)
        return full_source + translation
    rotation, translation = _kabsch_align(retained_source, retained_target)
    return (rotation @ full_source.T).T + translation


def _kabsch_align(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the rigid transform aligning source points to target points."""
    centroid_source = source.mean(axis=0)
    centroid_target = target.mean(axis=0)
    source_centered = source - centroid_source
    target_centered = target - centroid_target
    covariance = source_centered.T @ target_centered
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    determinant = np.linalg.det(vt_matrix.T @ u_matrix.T)
    correction = np.diag([1.0, 1.0, determinant])
    rotation = vt_matrix.T @ correction @ u_matrix.T
    translation = centroid_target - rotation @ centroid_source
    return rotation, translation


def _snap_reactive_atom_to_bond_length(
    coords: np.ndarray,
    reactive_atom: PdbAtomRecord,
    protein_link_atom: PdbAtomRecord,
    target_bond_length: float,
) -> np.ndarray:
    """Translate the modifier so the reactive atom sits at bond length."""
    if reactive_atom.atom_index is None:
        raise ValueError("Modifier reactive atom requires atom_index for coordinate snapping")
    site_coord = _coord(protein_link_atom)
    current_reactive = coords[reactive_atom.atom_index]
    direction = current_reactive - site_coord
    distance = np.linalg.norm(direction)
    if distance < 1e-8:
        raise RuntimeError("Modifier reactive atom coincides with the protein linker atom")
    target_reactive = site_coord + target_bond_length * direction / distance
    return coords + (target_reactive - current_reactive)


def _placed_fragment_from_coords(
    modifier: GeneratedPolymerFragment, coords: np.ndarray
) -> PlacedPolymerFragment:
    """Build a placed fragment preserving atom identity and connectivity."""
    placed_atoms = []
    for atom in modifier.to_placed_fragment().atoms:
        if atom.atom_index is None:
            raise ValueError("Modifier atom_index is required to restore placed coordinates")
        x_coord, y_coord, z_coord = coords[atom.atom_index]
        placed_atoms.append(
            atom.model_copy(update={"x": float(x_coord), "y": float(y_coord), "z": float(z_coord)})
        )
    placed = modifier.to_placed_fragment()
    return placed.model_copy(update={"atoms": tuple(placed_atoms)})


def _minimum_distance(points_a: np.ndarray, points_b: np.ndarray) -> float:
    """Return the minimum pairwise distance between two coordinate arrays."""
    if len(points_a) == 0 or len(points_b) == 0:
        return float("inf")
    deltas = points_a[:, np.newaxis, :] - points_b[np.newaxis, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    return float(np.min(distances))


def _same_atom(atom_a: PdbAtomRecord, atom_b: PdbAtomRecord) -> bool:
    """Return whether two PDB atom records have the same source identity."""
    if atom_a.atom_index is not None and atom_b.atom_index is not None:
        return atom_a.atom_index == atom_b.atom_index
    if atom_a.serial is not None and atom_b.serial is not None:
        return atom_a.serial == atom_b.serial
    return atom_a.atom_name.upper() == atom_b.atom_name.upper()


def _atom_identity(atom: PdbAtomRecord) -> tuple[int | None, int | None, str, int, str, str, str]:
    """Return a stable PDB atom identity for placement filtering."""
    return (
        atom.atom_index,
        atom.serial,
        atom.atom_name.upper(),
        atom.residue_number,
        atom.insertion_code.upper(),
        atom.residue_name.upper(),
        atom.chain_id.upper(),
    )
