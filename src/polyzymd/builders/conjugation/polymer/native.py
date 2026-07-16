"""Native bundled methacrylate polymer generation."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.builders.conjugation.polymer.mbuild import from_mbuild
from polyzymd.builders.conjugation.polymer.recipe import PolymerRecipe
from polyzymd.utils.charging import get_charger

NATIVE_METHACRYLATE_CACHE_VERSION = "native-methacrylate-v1"
INTER_RESIDUE_BOND_LENGTH_ANGSTROM = 1.54
HEAVY_CLASH_THRESHOLD_ANGSTROM = 1.0
MAX_PLACEMENT_RETRIES = 200


@dataclass(frozen=True)
class MethacrylateMatch:
    """Exact methacrylate alkene atom roles in an explicit-hydrogen monomer."""

    head_index: int
    tail_index: int


@dataclass(frozen=True)
class NativeGenerationPaths:
    """Generated artifact paths for one native methacrylate polymer."""

    pdb_path: Path
    sdf_path: Path
    charged_sdf_path: Path
    identity_path: Path


@dataclass(frozen=True)
class NativeMethacrylateResult:
    """In-memory and on-disk result from native methacrylate generation."""

    compound: Any
    molecule: Any
    charged_molecule: Any
    paths: NativeGenerationPaths
    sequence: str


def native_cache_directory(cache_directory: Path | str) -> Path:
    """Return the versioned native methacrylate cache directory.

    Parameters
    ----------
    cache_directory : pathlib.Path or str
        User-facing polymer cache root.

    Returns
    -------
    pathlib.Path
        Versioned subdirectory used only by the bundled native generator.
    """
    return Path(cache_directory) / NATIVE_METHACRYLATE_CACHE_VERSION


def native_cache_identity(
    recipe: PolymerRecipe, sequence: str, *, charger_type: str = "nagl"
) -> str:
    """Return a deterministic native methacrylate cache fingerprint.

    Parameters
    ----------
    recipe : PolymerRecipe
        Polymer recipe used for generation.
    sequence : str
        Exact sequence being generated.
    charger_type : str, optional
        Charge assignment method identity, by default ``"nagl"``.

    Returns
    -------
    str
        SHA256 hex digest covering chemistry, sequence, settings, and charge method.
    """
    payload = _cache_identity_payload(recipe, sequence, charger_type=charger_type)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def native_artifact_paths(
    recipe: PolymerRecipe,
    sequence: str,
    cache_directory: Path | str,
    *,
    charger_type: str = "nagl",
) -> NativeGenerationPaths:
    """Return centralized native artifact paths for a generated sequence.

    Parameters
    ----------
    recipe : PolymerRecipe
        Recipe defining monomer identities.
    sequence : str
        Exact sequence labels.
    cache_directory : pathlib.Path or str
        User-facing polymer cache root.
    charger_type : str, optional
        Charge assignment method, by default ``"nagl"``.

    Returns
    -------
    NativeGenerationPaths
        Versioned and fingerprinted cache artifact paths.
    """
    identity = native_cache_identity(recipe, sequence, charger_type=charger_type)
    directory = native_cache_directory(cache_directory) / identity[:16]
    stem = _artifact_stem(recipe, sequence)
    return NativeGenerationPaths(
        pdb_path=directory / f"{stem}.pdb",
        sdf_path=directory / f"{stem}.sdf",
        charged_sdf_path=directory / f"{stem}_charged.sdf",
        identity_path=directory / "native_identity.json",
    )


def methacrylate_match(smiles: str) -> MethacrylateMatch:
    """Find the single exact methacrylate alkene in a monomer SMILES string.

    The bundled recipe is intentionally narrow. It accepts exactly one C=C where
    one carbon is terminal CH2 and the substituted carbon has no hydrogens, one
    carbonyl-carbon neighbor, and one methyl-carbon neighbor.

    Parameters
    ----------
    smiles : str
        Monomer SMILES string.

    Returns
    -------
    MethacrylateMatch
        Head/tail atom indices in an explicit-hydrogen RDKit molecule.

    Raises
    ------
    ValueError
        If the monomer has zero or multiple exact methacrylate matches.
    """
    from rdkit import Chem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        raise ValueError("Could not parse monomer SMILES for native methacrylate generation")
    head, tail = _methacrylate_match_in_mol(mol)
    return MethacrylateMatch(head_index=head, tail_index=tail)


def generate_native_methacrylate_polymer(
    recipe: PolymerRecipe,
    cache_directory: Path | str,
    *,
    sequence: str | None = None,
    charger_type: str = "nagl",
    force_regenerate: bool = False,
) -> NativeMethacrylateResult:
    """Generate a polymer with the bundled native methacrylate recipe.

    Parameters
    ----------
    recipe : PolymerRecipe
        Validated recipe containing methacrylate monomer SMILES.
    cache_directory : pathlib.Path or str
        User-facing polymer cache root. Native artifacts are written under a
        versioned subdirectory to avoid reusing older backend outputs.
    sequence : str or None, optional
        Exact sequence to build. When omitted, ``recipe.generate_sequence()`` is
        used, by default ``None``.
    charger_type : str, optional
        Charge assignment backend resolved through ``get_charger()``, by default
        ``"nagl"``.
    force_regenerate : bool, optional
        Rebuild even when the charged native SDF exists, by default ``False``.

    Returns
    -------
    NativeMethacrylateResult
        Generated mBuild compound, OpenFF molecules, paths, and sequence.
    """
    from openff.toolkit import Molecule

    selected_sequence = sequence or recipe.generate_sequence()
    paths = native_artifact_paths(
        recipe,
        selected_sequence,
        cache_directory,
        charger_type=charger_type,
    )
    paths.charged_sdf_path.parent.mkdir(parents=True, exist_ok=True)

    if paths.charged_sdf_path.exists() and not force_regenerate:
        cached_charged = Molecule.from_file(
            paths.charged_sdf_path,
            file_format="SDF",
            allow_undefined_stereo=True,
        )
        compound = assemble_native_methacrylate_compound(recipe, selected_sequence)
        molecule = from_mbuild(compound)
        if cached_charged.n_atoms != molecule.n_atoms:
            raise ValueError("Cached native methacrylate SDF atom count does not match the recipe")
        molecule.partial_charges = cached_charged.partial_charges
        return NativeMethacrylateResult(
            compound=compound,
            molecule=molecule,
            charged_molecule=molecule,
            paths=paths,
            sequence=selected_sequence,
        )

    compound = assemble_native_methacrylate_compound(recipe, selected_sequence)
    molecule = from_mbuild(compound)
    _validate_openff_geometry(molecule)
    _atomic_write_json(
        _cache_identity_payload(recipe, selected_sequence, charger_type=charger_type),
        paths.identity_path,
    )
    _atomic_write_openff_pdb(molecule, paths.pdb_path)
    _atomic_write_openff_file(molecule, paths.sdf_path)

    charged = get_charger(charger_type).charge_molecule(molecule)
    _validate_openff_geometry(charged)
    _atomic_write_openff_file(charged, paths.charged_sdf_path)

    return NativeMethacrylateResult(
        compound=compound,
        molecule=molecule,
        charged_molecule=charged,
        paths=paths,
        sequence=selected_sequence,
    )


def assemble_native_methacrylate_compound(recipe: PolymerRecipe, sequence: str) -> Any:
    """Assemble a native mBuild compound for a recipe-defined sequence.

    Parameters
    ----------
    recipe : PolymerRecipe
        Recipe supplying monomer SMILES, names, and residue names.
    sequence : str
        Exact sequence labels to assemble.

    Returns
    -------
    mbuild.Compound
        Hierarchical compound with one residue child per sequence position.
    """
    if len(sequence) != recipe.length:
        raise ValueError("Native methacrylate sequence length must match the recipe length")
    if len(sequence) < 2:
        raise ValueError("Native methacrylate generation requires at least two monomers")

    from rdkit import Chem

    monomers: list[tuple[str, Any, int, int]] = []
    for residue_index, label in enumerate(sequence):
        monomer = recipe.monomer_by_label(label)
        if residue_index in {0, len(sequence) - 1}:
            monomers.append(_terminal_variant(monomer.residue_name, monomer.smiles))
        else:
            monomers.append(_internal_variant(monomer.residue_name, monomer.smiles))

    combo = Chem.RWMol()
    offsets: list[int] = []
    alkene_pairs: list[tuple[int, int]] = []
    conformer_positions: list[np.ndarray] = []
    for _residue_name, mol, head, tail in monomers:
        offsets.append(combo.GetNumAtoms())
        combo.InsertMol(mol)
        alkene_pairs.append((offsets[-1] + head, offsets[-1] + tail))
        conformer_positions.append(_conformer_positions(mol))

    combo.AddBond(alkene_pairs[0][0], alkene_pairs[1][0], Chem.BondType.SINGLE)
    for index in range(1, len(sequence) - 1):
        combo.AddBond(alkene_pairs[index][1], alkene_pairs[index + 1][0], Chem.BondType.SINGLE)

    mol = combo.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    _set_random_walk_conformer(
        mol, monomers, offsets, alkene_pairs, conformer_positions, recipe, sequence
    )
    return _mbuild_compound_from_rdkit(mol, monomers, offsets, sequence)


def _terminal_variant(residue_name: str, smiles: str) -> tuple[str, Any, int, int]:
    """Return a one-port terminal methacrylate variant."""
    from rdkit import Chem

    editable = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(smiles)))
    head, tail = _methacrylate_match_in_mol(editable.GetMol())
    editable.RemoveAtom(_one_hydrogen_neighbor_index(editable, head))
    mol = editable.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    _embed_variant_conformer(mol, f"terminal:{residue_name}:{smiles}")
    return residue_name, mol, head, tail


def _internal_variant(residue_name: str, smiles: str) -> tuple[str, Any, int, int]:
    """Return a two-port saturated internal methacrylate variant."""
    from rdkit import Chem

    editable = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(smiles)))
    head, tail = _methacrylate_match_in_mol(editable.GetMol())
    editable.RemoveBond(head, tail)
    editable.AddBond(head, tail, Chem.BondType.SINGLE)
    mol = editable.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    _embed_variant_conformer(mol, f"internal:{residue_name}:{smiles}")
    return residue_name, mol, head, tail


def _methacrylate_match_in_mol(mol: Any) -> tuple[int, int]:
    """Return the unique native methacrylate head and tail indices."""
    matches = []
    for bond in mol.GetBonds():
        if float(bond.GetBondTypeAsDouble()) != 2.0:
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetSymbol() != "C" or end.GetSymbol() != "C":
            continue
        ordered = _ordered_methacrylate_pair(begin, end)
        if ordered is not None:
            matches.append(ordered)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one methacrylate alkene match, found {len(matches)}")
    return matches[0]


def _ordered_methacrylate_pair(atom_a: Any, atom_b: Any) -> tuple[int, int] | None:
    """Return head/tail indices when two alkene atoms match the bundled pattern."""
    for head, tail in ((atom_a, atom_b), (atom_b, atom_a)):
        if _hydrogen_neighbor_count(head) != 2 or _hydrogen_neighbor_count(tail) != 0:
            continue
        carbon_neighbors = [
            neighbor
            for neighbor in tail.GetNeighbors()
            if neighbor.GetIdx() != head.GetIdx() and neighbor.GetSymbol() == "C"
        ]
        if len(carbon_neighbors) != 2:
            continue
        if not any(_is_carbonyl_carbon(neighbor, tail.GetIdx()) for neighbor in carbon_neighbors):
            continue
        if not any(_hydrogen_neighbor_count(neighbor) == 3 for neighbor in carbon_neighbors):
            continue
        return head.GetIdx(), tail.GetIdx()
    return None


def _hydrogen_neighbor_count(atom: Any) -> int:
    """Return the number of explicit hydrogen neighbors on an atom."""
    return sum(neighbor.GetSymbol() == "H" for neighbor in atom.GetNeighbors())


def _is_carbonyl_carbon(atom: Any, excluded_neighbor_index: int) -> bool:
    """Return whether an atom has a double-bonded oxygen neighbor."""
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetIdx() == excluded_neighbor_index:
            continue
        if other.GetSymbol() == "O" and float(bond.GetBondTypeAsDouble()) == 2.0:
            return True
    return False


def _one_hydrogen_neighbor_index(editable_mol: Any, atom_index: int) -> int:
    """Return one removable explicit hydrogen neighbor index."""
    atom = editable_mol.GetAtomWithIdx(atom_index)
    hydrogen_indices = sorted(
        neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "H"
    )
    if not hydrogen_indices:
        raise ValueError("Terminal methacrylate vinyl atom has no removable hydrogen")
    return hydrogen_indices[-1]


def _embed_variant_conformer(mol: Any, identity: str) -> None:
    """Generate a deterministic chemically reasonable conformer for a variant."""
    from rdkit.Chem import AllChem

    seed = int(hashlib.sha256(identity.encode("utf-8")).hexdigest()[:8], 16) % (2**31 - 1)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = False
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise ValueError(f"RDKit could not embed native methacrylate variant {identity!r}")
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except (RuntimeError, ValueError):
        # Coordinates from ETKDG are kept when UFF lacks parameters
        return


def _conformer_positions(mol: Any) -> np.ndarray:
    """Return conformer positions from an RDKit molecule as an array."""
    conformer = mol.GetConformer()
    return np.array(
        [
            [
                conformer.GetAtomPosition(index).x,
                conformer.GetAtomPosition(index).y,
                conformer.GetAtomPosition(index).z,
            ]
            for index in range(mol.GetNumAtoms())
        ],
        dtype=float,
    )


def _set_random_walk_conformer(
    mol: Any,
    monomers: list[tuple[str, Any, int, int]],
    offsets: list[int],
    alkene_pairs: list[tuple[int, int]],
    conformer_positions: list[np.ndarray],
    recipe: PolymerRecipe,
    sequence: str,
) -> None:
    """Place repeat-unit conformers with a deterministic self-avoiding random walk."""
    from rdkit import Chem

    rng = np.random.default_rng(int(native_cache_identity(recipe, sequence)[:16], 16))
    placed = np.full((mol.GetNumAtoms(), 3), np.nan, dtype=float)
    placed_heavy_indices: list[int] = []

    for residue_index, ((_residue_name, residue_mol, head, tail), offset, positions) in enumerate(
        zip(monomers, offsets, conformer_positions, strict=True)
    ):
        atom_indices = list(range(offset, offset + residue_mol.GetNumAtoms()))
        if residue_index == 0:
            transformed = positions - positions[head]
        else:
            previous_port = alkene_pairs[residue_index - 1][0]
            if residue_index > 1:
                previous_port = alkene_pairs[residue_index - 1][1]
            previous_position = placed[previous_port]
            transformed = _place_residue_with_retries(
                positions,
                residue_mol,
                head,
                previous_position,
                placed,
                placed_heavy_indices,
                rng,
            )
        for local_index, global_index in enumerate(atom_indices):
            placed[global_index] = transformed[local_index]
            if residue_mol.GetAtomWithIdx(local_index).GetSymbol() != "H":
                placed_heavy_indices.append(global_index)

    if not np.isfinite(placed).all():
        raise ValueError(
            "Native methacrylate random-walk placement produced non-finite coordinates"
        )

    conformer = Chem.Conformer(mol.GetNumAtoms())
    for index, (x_coord, y_coord, z_coord) in enumerate(placed):
        conformer.SetAtomPosition(index, (float(x_coord), float(y_coord), float(z_coord)))
    mol.RemoveAllConformers()
    mol.AddConformer(conformer, assignId=True)


def _place_residue_with_retries(
    positions: np.ndarray,
    residue_mol: Any,
    head_index: int,
    previous_position: np.ndarray,
    placed: np.ndarray,
    placed_heavy_indices: list[int],
    rng: np.random.Generator,
) -> np.ndarray:
    """Place one residue by random orientation with heavy-atom clash rejection."""
    centered = positions - positions[head_index]
    best_candidate = None
    best_clash = -np.inf
    for _attempt in range(MAX_PLACEMENT_RETRIES):
        direction = _random_unit_vector(rng)
        target = previous_position + direction * INTER_RESIDUE_BOND_LENGTH_ANGSTROM
        candidate = centered @ _random_rotation_matrix(rng).T + target
        clash_distance = _minimum_new_heavy_distance(
            candidate,
            residue_mol,
            placed,
            placed_heavy_indices,
        )
        if clash_distance >= HEAVY_CLASH_THRESHOLD_ANGSTROM:
            return candidate
        if clash_distance > best_clash:
            best_clash = clash_distance
            best_candidate = candidate
    raise ValueError(
        "Native methacrylate random-walk placement failed after "
        f"{MAX_PLACEMENT_RETRIES} retries; closest accepted heavy distance would be "
        f"{best_clash:.2f} Å, below the {HEAVY_CLASH_THRESHOLD_ANGSTROM:.2f} Å threshold"
    ) from None


def _minimum_new_heavy_distance(
    candidate: np.ndarray,
    residue_mol: Any,
    placed: np.ndarray,
    placed_heavy_indices: list[int],
) -> float:
    """Return the minimum distance from new heavy atoms to already placed heavy atoms."""
    if not placed_heavy_indices:
        return float("inf")
    new_heavy = np.array(
        [
            candidate[index]
            for index, atom in enumerate(residue_mol.GetAtoms())
            if atom.GetSymbol() != "H"
        ]
    )
    old_heavy = placed[placed_heavy_indices]
    distances = np.linalg.norm(new_heavy[:, None, :] - old_heavy[None, :, :], axis=2)
    return float(np.min(distances))


def _random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    """Return one deterministic random unit vector."""
    vector = rng.normal(size=3)
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        return np.array([1.0, 0.0, 0.0])
    return vector / norm


def _random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """Return a uniformly sampled deterministic 3D rotation matrix."""
    q1, q2, q3 = rng.random(3)
    quat = np.array(
        [
            np.sqrt(1 - q1) * np.sin(2 * np.pi * q2),
            np.sqrt(1 - q1) * np.cos(2 * np.pi * q2),
            np.sqrt(q1) * np.sin(2 * np.pi * q3),
            np.sqrt(q1) * np.cos(2 * np.pi * q3),
        ]
    )
    x_val, y_val, z_val, w_val = quat
    return np.array(
        [
            [
                1 - 2 * (y_val**2 + z_val**2),
                2 * (x_val * y_val - z_val * w_val),
                2 * (x_val * z_val + y_val * w_val),
            ],
            [
                2 * (x_val * y_val + z_val * w_val),
                1 - 2 * (x_val**2 + z_val**2),
                2 * (y_val * z_val - x_val * w_val),
            ],
            [
                2 * (x_val * z_val - y_val * w_val),
                2 * (y_val * z_val + x_val * w_val),
                1 - 2 * (x_val**2 + y_val**2),
            ],
        ]
    )


def _mbuild_compound_from_rdkit(
    mol: Any,
    monomers: list[tuple[str, Any, int, int]],
    offsets: list[int],
    sequence: str,
) -> Any:
    """Build an mBuild hierarchy from an RDKit molecule and residue spans."""
    import mbuild as mb

    compound = mb.Compound(name="native-methacrylate")
    particles = []
    residue_by_atom = _residue_by_atom(monomers, offsets, mol.GetNumAtoms(), sequence)
    residue_compounds = {}
    for atom in mol.GetAtoms():
        residue_name, residue_number = residue_by_atom[atom.GetIdx()]
        key = (residue_name, residue_number)
        if key not in residue_compounds:
            residue = mb.Compound(name=residue_name)
            residue.residue_number = residue_number
            compound.add(residue)
            residue_compounds[key] = residue
        coord = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        particle = mb.Compound(
            name=f"{atom.GetSymbol()}{atom.GetIdx() + 1}",
            pos=np.array([coord.x, coord.y, coord.z], dtype=float) / 10.0,
        )
        particle.element = atom.GetSymbol()
        particle.formal_charge = atom.GetFormalCharge()
        residue_compounds[key].add(particle)
        particles.append(particle)
    for bond in mol.GetBonds():
        compound.add_bond(
            (particles[bond.GetBeginAtomIdx()], particles[bond.GetEndAtomIdx()]),
            bond_order=float(bond.GetBondTypeAsDouble()),
        )
    return compound


def _residue_by_atom(
    monomers: list[tuple[str, Any, int, int]],
    offsets: list[int],
    atom_count: int,
    sequence: str,
) -> dict[int, tuple[str, int]]:
    """Return residue metadata for assembled monomer atoms."""
    mapping = {}
    for residue_number, ((residue_name, mol, _head, _tail), offset) in enumerate(
        zip(monomers, offsets, strict=True),
        start=1,
    ):
        for atom in mol.GetAtoms():
            mapping[offset + atom.GetIdx()] = (residue_name, residue_number)
    for atom_index in range(atom_count):
        mapping.setdefault(atom_index, ("CAP", len(sequence) + 1))
    return mapping


def _artifact_stem(recipe: PolymerRecipe, sequence: str) -> str:
    """Return the native artifact stem matching the legacy filename shape."""
    monomer_names = recipe.to_sequence_monomer_names()
    unique_labels = []
    for label in sequence:
        if label not in unique_labels:
            unique_labels.append(label)
    monomer_prefix = "-".join(monomer_names[label] for label in unique_labels)
    return f"{monomer_prefix}_seq={sequence}_{len(sequence)}-mer"


def _cache_identity_payload(
    recipe: PolymerRecipe,
    sequence: str,
    *,
    charger_type: str,
) -> dict[str, Any]:
    """Return the canonical native cache identity payload."""
    return {
        "mechanism": "bundled-native-methacrylate",
        "schema_version": NATIVE_METHACRYLATE_CACHE_VERSION,
        "sequence": sequence,
        "recipe": {
            "name": recipe.name,
            "length": recipe.length,
            "seed": recipe.seed,
            "reactive_monomer_label": recipe.reactive_monomer_label,
            "reactive_monomer_index": recipe.reactive_monomer_index,
            "fixed_sequence": recipe.fixed_sequence,
            "monomers": [
                {
                    "label": monomer.label,
                    "name": monomer.name,
                    "residue_name": monomer.residue_name,
                    "probability": monomer.probability,
                    "smiles": monomer.smiles,
                    "canonical_smiles": _canonical_smiles(monomer.smiles),
                }
                for monomer in recipe.monomers
            ],
        },
        "charge": _charge_identity(charger_type),
    }


def _canonical_smiles(smiles: str) -> str:
    """Return an RDKit canonical SMILES string for cache identity."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Could not parse monomer SMILES for cache identity")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _charge_identity(charger_type: str) -> dict[str, str]:
    """Return cache identity fields for the selected charge method."""
    normalized = charger_type.lower().strip()
    if normalized == "nagl":
        from polyzymd.utils.charging import DEFAULT_NAGL_MODEL

        return {"method": normalized, "model": DEFAULT_NAGL_MODEL}
    return {"method": normalized, "model": "default"}


def _atomic_write_json(payload: dict[str, Any], path: Path) -> None:
    """Atomically write a JSON payload in the destination directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, delete=False, encoding="utf-8"
    ) as handle:
        temp_path = Path(handle.name)
        json.dump(payload, handle, sort_keys=True, indent=2)
        handle.write("\n")
    os.replace(temp_path, path)


def _atomic_write_openff_file(molecule: Any, path: Path) -> None:
    """Atomically write an OpenFF molecule file in the destination directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=path.suffix, delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        molecule.to_file(temp_path, file_format="SDF")
        os.replace(temp_path, path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _atomic_write_openff_pdb(molecule: Any, path: Path) -> None:
    """Atomically write a minimal PDB with residue metadata and CONECT records."""
    coordinates = molecule.conformers[0].m_as("angstrom")
    lines = []
    for index, atom in enumerate(molecule.atoms, start=1):
        metadata = atom.metadata
        atom_name = str(metadata.get("atom_name", f"{atom.symbol}{index}"))[:4]
        residue_name = str(metadata.get("residue_name", "MOL"))[:3]
        residue_number = int(metadata.get("residue_number", 1))
        x, y, z = coordinates[index - 1]
        lines.append(
            f"HETATM{index:5d} {atom_name:<4} {residue_name:>3} C{residue_number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom.symbol:>2}\n"
        )
    for bond in molecule.bonds:
        lines.append(f"CONECT{bond.atom1_index + 1:5d}{bond.atom2_index + 1:5d}\n")
    lines.append("END\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, delete=False, encoding="utf-8"
    ) as handle:
        temp_path = Path(handle.name)
        handle.write("".join(lines))
    os.replace(temp_path, path)


def _validate_openff_geometry(molecule: Any) -> None:
    """Validate finite, non-collinear, and clash-free OpenFF coordinates."""
    coords = molecule.conformers[0].m_as("angstrom")
    if not np.isfinite(coords).all():
        raise ValueError("Native methacrylate coordinates contain non-finite values")
    if np.linalg.matrix_rank(coords - np.mean(coords, axis=0), tol=1.0e-3) < 2:
        raise ValueError("Native methacrylate coordinates are collinear or degenerate")
    for bond in molecule.bonds:
        distance = float(np.linalg.norm(coords[bond.atom1_index] - coords[bond.atom2_index]))
        if not 0.7 <= distance <= 2.2:
            raise ValueError(f"Native methacrylate bond length {distance:.2f} Å is implausible")
    graph_distances = _graph_distances_within_two_bonds(molecule)
    heavy_indices = [atom.molecule_atom_index for atom in molecule.atoms if atom.symbol != "H"]
    for position, left in enumerate(heavy_indices):
        for right in heavy_indices[position + 1 :]:
            if tuple(sorted((left, right))) in graph_distances:
                continue
            distance = float(np.linalg.norm(coords[left] - coords[right]))
            if distance < HEAVY_CLASH_THRESHOLD_ANGSTROM:
                raise ValueError(
                    "Native methacrylate heavy-atom clash detected: "
                    f"atoms {left + 1}-{right + 1} are {distance:.2f} Å apart"
                )


def _graph_distances_within_two_bonds(molecule: Any) -> set[tuple[int, int]]:
    """Return atom pairs separated by at most two graph bonds."""
    neighbors: dict[int, set[int]] = {atom.molecule_atom_index: set() for atom in molecule.atoms}
    for bond in molecule.bonds:
        neighbors[bond.atom1_index].add(bond.atom2_index)
        neighbors[bond.atom2_index].add(bond.atom1_index)
    close_pairs = set()
    for atom_index, bonded in neighbors.items():
        for other in bonded:
            close_pairs.add(tuple(sorted((atom_index, other))))
            for second in neighbors[other]:
                if second != atom_index:
                    close_pairs.add(tuple(sorted((atom_index, second))))
    return close_pairs
