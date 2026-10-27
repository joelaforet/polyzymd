"""Prepare a narrow 4CHA PDB example for OpenFF ingestion checks.

This script selects the first 4CHA alpha-chymotrypsin enzyme copy, removes
heterogens and waters, adds missing hydrogens only, and relabels retained protein
chains to PolyzyMD chain ``A``. It does not model missing residues or missing
heavy atoms.
"""

from __future__ import annotations

import argparse
from pathlib import Path

KEEP_CHAINS = {"A", "B", "C"}


def prepare_4cha(
    raw_pdb: Path,
    output_pdb: Path,
    ph: float,
    rename_first_cys_to_ncyx: bool,
) -> None:
    """Prepare one 4CHA enzyme copy.

    Heavy simulation-stack imports stay inside the function so importing this
    script does not require OpenMM or PDBFixer.
    """
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(raw_pdb))

    chains_to_remove = [
        chain_index
        for chain_index, chain in enumerate(fixer.topology.chains())
        if chain.id not in KEEP_CHAINS
    ]
    fixer.removeChains(chains_to_remove)
    fixer.removeHeterogens(keepWater=False)

    # Deliberately add hydrogens only. Missing residues/heavy atoms are external
    # modeling decisions and should be reviewed before simulation.
    fixer.addMissingHydrogens(ph)

    if rename_first_cys_to_ncyx:
        first_residue = next(fixer.topology.residues(), None)
        if first_residue is None or first_residue.name != "CYS":
            raise ValueError("Expected the first retained 4CHA residue to be CYS")
        first_residue.name = "NCYX"
        print("Renamed first retained CYS residue to NCYX for the proof of concept")

    for chain in fixer.topology.chains():
        chain.id = "A"

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    print(f"Wrote {output_pdb}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_pdb", type=Path, help="Raw 4CHA PDB path")
    parser.add_argument("output_pdb", type=Path, help="Prepared output PDB path")
    parser.add_argument("--ph", type=float, default=7.0, help="pH for PDBFixer hydrogens")
    parser.add_argument(
        "--rename-first-cys-to-ncyx",
        action="store_true",
        help="Rename the first retained CYS residue to NCYX for the private-API proof of concept",
    )
    return parser.parse_args()


def main() -> None:
    """Run the example preparation."""
    args = parse_args()
    prepare_4cha(args.raw_pdb, args.output_pdb, args.ph, args.rename_first_cys_to_ncyx)


if __name__ == "__main__":
    main()
