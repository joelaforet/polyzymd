"""Protein preparation helpers for OpenFF/Pablo conjugation workflows."""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import BaseModel, Field


class ProteinCanonicalizationSettings(BaseModel):
    """Settings for source-protein hydrogen canonicalization."""

    ph: float = Field(7.0, description="pH used when OpenMM regenerates protein hydrogens")
    force_field_name: str = Field(
        "amber14/protein.ff14SB.xml",
        description="OpenMM force field XML used by Modeller.addHydrogens",
    )
    output_path: Path | None = Field(
        None,
        description="Optional canonicalized PDB output path",
    )


class ProteinCanonicalizationResult(BaseModel):
    """Result of source-protein hydrogen canonicalization."""

    input_path: Path
    output_path: Path
    ph: float
    force_field_name: str
    atom_count: int
    hydrogen_count: int
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    def save(self, path: Path | str | None = None) -> Path:
        """Write a JSON summary and return the output path."""
        output_path = Path(path) if path is not None else self.output_path.with_suffix(".json")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")
        return output_path


def canonicalize_protein_hydrogens(
    pdb_path: Path | str,
    output_path: Path | str | None = None,
    *,
    settings: ProteinCanonicalizationSettings | None = None,
) -> ProteinCanonicalizationResult:
    """Regenerate source-protein hydrogens with PDBFixer and OpenMM Modeller.

    PDBFixer alone can preserve noncanonical hydrogen names from input PDBs.
    This helper intentionally deletes all existing hydrogens after PDBFixer
    preprocessing and asks OpenMM Modeller to add pH-specific canonical protein
    hydrogens. Chain IDs and residue numbers are preserved when OpenMM can keep
    the original IDs.
    """
    preparation_settings = settings or ProteinCanonicalizationSettings()
    input_path = Path(pdb_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Protein PDB not found: {input_path}")

    destination = Path(output_path) if output_path is not None else preparation_settings.output_path
    if destination is None:
        destination = input_path.with_name(
            f"{input_path.stem}_canonical_h_pH{preparation_settings.ph:g}.pdb"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)

    try:
        from openmm.app import ForceField, Modeller, PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "Protein hydrogen canonicalization requires PDBFixer and OpenMM; "
            "run in the PolyzyMD build pixi environment."
        ) from exc

    fixer = PDBFixer(filename=str(input_path))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    modeller = Modeller(fixer.topology, fixer.positions)
    hydrogens = [atom for atom in modeller.topology.atoms() if atom.element.symbol == "H"]
    if hydrogens:
        modeller.delete(hydrogens)

    force_field = ForceField(preparation_settings.force_field_name)
    modeller.addHydrogens(force_field, pH=preparation_settings.ph)

    with destination.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)

    atom_count, hydrogen_count = _count_pdb_atoms(destination)
    return ProteinCanonicalizationResult(
        input_path=input_path,
        output_path=destination,
        ph=preparation_settings.ph,
        force_field_name=preparation_settings.force_field_name,
        atom_count=atom_count,
        hydrogen_count=hydrogen_count,
        diagnostics=(
            "Source protein hydrogens were removed and regenerated with OpenMM Modeller",
            f"Hydrogen canonicalization pH: {preparation_settings.ph:g}",
            f"OpenMM hydrogen force field: {preparation_settings.force_field_name}",
        ),
    )


def _count_pdb_atoms(path: Path) -> tuple[int, int]:
    """Return total atom and hydrogen counts for a PDB file."""
    atom_count = 0
    hydrogen_count = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_count += 1
        element = line[76:78].strip()
        atom_name = line[12:16].strip()
        if element.upper() == "H" or atom_name.upper().startswith("H"):
            hydrogen_count += 1
    return atom_count, hydrogen_count


__all__ = [
    "ProteinCanonicalizationResult",
    "ProteinCanonicalizationSettings",
    "canonicalize_protein_hydrogens",
]
