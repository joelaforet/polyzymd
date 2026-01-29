#!/usr/bin/env python
"""Generator script for pre-computed solvent SDF files.

This script generates SDF files with embedded AM1BCC partial charges for
all solvents in the cosolvent library. These files are shipped with the
package to provide instant loading of common solvents.

Usage:
    python -m polyzymd.data.solvents._generator

    Or from the package root:
    python src/polyzymd/data/solvents/_generator.py

The script will:
1. Generate 3D conformers for each solvent
2. Compute AM1BCC partial charges
3. Save to SDF files with embedded charges and metadata

This only needs to be run once during development. The resulting SDF files
are committed to the repository and shipped with the package.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
LOGGER = logging.getLogger(__name__)


def generate_all_solvents() -> None:
    """Generate SDF files for all library solvents."""
    from openff.toolkit import Molecule
    from openff.units import unit as offunit

    from polyzymd.data.cosolvent_library import COSOLVENT_LIBRARY

    # Output directory (same as this script's directory)
    output_dir = Path(__file__).parent
    LOGGER.info(f"Output directory: {output_dir}")

    # Track results
    successful = []
    failed = []

    # Generate TIP3P water first
    LOGGER.info("=" * 60)
    LOGGER.info("Generating TIP3P water...")
    try:
        water = _generate_tip3p_water()
        output_path = output_dir / "tip3p.sdf"
        _save_molecule(water, output_path)
        successful.append("tip3p")
        LOGGER.info(f"  ✓ Saved to {output_path}")
    except Exception as e:
        LOGGER.error(f"  ✗ Failed: {e}")
        failed.append(("tip3p", str(e)))

    # Generate each cosolvent
    for name, data in COSOLVENT_LIBRARY.items():
        LOGGER.info("=" * 60)
        LOGGER.info(f"Generating {name} ({data.name})...")
        LOGGER.info(f"  SMILES: {data.smiles}")

        try:
            mol = _generate_charged_molecule(
                smiles=data.smiles,
                residue_name=name[:3].upper(),
                name=name,
            )

            output_path = output_dir / f"{name}.sdf"
            _save_molecule(mol, output_path)

            # Report charges
            charges = mol.partial_charges.magnitude
            LOGGER.info(f"  Total charge: {sum(charges):.6f}")
            LOGGER.info(f"  Charge range: [{min(charges):.4f}, {max(charges):.4f}]")
            LOGGER.info(f"  ✓ Saved to {output_path}")

            successful.append(name)

        except Exception as e:
            LOGGER.error(f"  ✗ Failed: {e}")
            failed.append((name, str(e)))

    # Summary
    LOGGER.info("=" * 60)
    LOGGER.info("SUMMARY")
    LOGGER.info("=" * 60)
    LOGGER.info(f"Successful: {len(successful)}")
    for name in successful:
        LOGGER.info(f"  ✓ {name}")

    if failed:
        LOGGER.info(f"Failed: {len(failed)}")
        for name, error in failed:
            LOGGER.info(f"  ✗ {name}: {error}")

    LOGGER.info("=" * 60)
    LOGGER.info("Done!")


def _generate_tip3p_water() -> "Molecule":
    """Generate TIP3P water with literature charges."""
    from openff.toolkit import Molecule
    from openff.units import unit as offunit

    TIP3P_CHARGES = {"O": -0.834, "H": 0.417}

    water = Molecule.from_smiles("O")
    water.name = "tip3p"

    # Generate 3D coordinates
    water.generate_conformers(n_conformers=1)

    # Assign literature charges
    charges = [TIP3P_CHARGES[atom.symbol] for atom in water.atoms]
    water.partial_charges = charges * offunit.elementary_charge

    # Set residue metadata
    for atom in water.atoms:
        atom.metadata["residue_name"] = "HOH"

    return water


def _generate_charged_molecule(
    smiles: str,
    residue_name: str,
    name: str,
) -> "Molecule":
    """Generate a molecule with AM1BCC charges."""
    from openff.toolkit import Molecule

    mol = Molecule.from_smiles(smiles)
    mol.name = name

    # Generate 3D conformer
    mol.generate_conformers(n_conformers=1)

    if mol.n_conformers == 0:
        raise RuntimeError(f"Failed to generate conformer for {name}")

    # Compute AM1BCC charges
    mol.assign_partial_charges(partial_charge_method="am1bcc")

    # Set residue metadata
    for atom in mol.atoms:
        atom.metadata["residue_name"] = residue_name

    return mol


def _save_molecule(mol: "Molecule", path: Path) -> None:
    """Save molecule to SDF with metadata."""
    import json

    # Package metadata
    metadata_dict = {}
    for i, atom in enumerate(mol.atoms):
        if atom.metadata:
            metadata_dict[str(i)] = dict(atom.metadata)

    if metadata_dict:
        mol.properties["metadata"] = json.dumps(metadata_dict)

    # Save to SDF
    mol.to_file(str(path), file_format="sdf")


def verify_sdf_files() -> None:
    """Verify that generated SDF files can be loaded correctly."""
    from openff.toolkit import Molecule

    solvents_dir = Path(__file__).parent

    LOGGER.info("Verifying generated SDF files...")

    for sdf_file in sorted(solvents_dir.glob("*.sdf")):
        try:
            mol = Molecule.from_file(str(sdf_file), file_format="sdf")

            # Check charges are present
            if mol.partial_charges is None:
                LOGGER.warning(f"  ⚠ {sdf_file.name}: No partial charges!")
            else:
                total_charge = sum(mol.partial_charges.magnitude)
                LOGGER.info(
                    f"  ✓ {sdf_file.name}: {mol.n_atoms} atoms, total charge = {total_charge:.6f}"
                )

        except Exception as e:
            LOGGER.error(f"  ✗ {sdf_file.name}: {e}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--verify":
        verify_sdf_files()
    else:
        generate_all_solvents()
        print()
        verify_sdf_files()
