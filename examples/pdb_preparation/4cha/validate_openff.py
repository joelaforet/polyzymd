"""Validate a prepared protein PDB with direct OpenFF ingestion.

This script is a diagnostic helper. It does not make a failed PDB valid; it
prints structural checks and then calls ``Topology.from_pdb`` directly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def summarize_pdb(pdb_path: Path) -> None:
    """Print lightweight structural checks before OpenFF validation."""
    lines = pdb_path.read_text(encoding="utf-8").splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    chains = sorted({line[21] for line in atoms})
    residues = sorted({line[17:20].strip() for line in atoms})
    elements = {line[76:78].strip() for line in atoms if len(line) >= 78}

    print("Structural checks")
    print(f"  atom records: {len(atoms)}")
    print(f"  chains: {chains}")
    print(f"  TER records: {sum(line.startswith('TER') for line in lines)}")
    print(f"  SSBOND records: {sum(line.startswith('SSBOND') for line in lines)}")
    print(f"  CONECT records: {sum(line.startswith('CONECT') for line in lines)}")
    print(f"  hydrogens present: {'H' in elements}")
    print(f"  residue names include CYX: {'CYX' in residues}")


def load_custom_substructures(path: Path | None) -> dict[str, dict[str, list[str]]] | None:
    """Load an OpenFF private custom-substructures mapping, if provided."""
    if path is None:
        return None
    with path.open(encoding="utf-8") as handle:
        data: Any = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("custom substructures JSON must be an object")
    return data


def validate_openff(
    pdb_path: Path,
    custom_substructures_path: Path | None = None,
) -> None:
    """Run direct OpenFF PDB ingestion validation."""
    from openff.toolkit import Topology

    summarize_pdb(pdb_path)
    custom_substructures = load_custom_substructures(custom_substructures_path)

    print("Running OpenFF Topology.from_pdb validation")
    if custom_substructures is None:
        topology = Topology.from_pdb(str(pdb_path))
    else:
        print(
            "  using private _custom_substructures proof of concept from "
            f"{custom_substructures_path}"
        )
        topology = Topology.from_pdb(
            str(pdb_path),
            _custom_substructures=custom_substructures,
        )

    print("OpenFF ingestion succeeded")
    print(f"  molecules: {topology.n_molecules}")
    print(f"  atoms: {topology.n_atoms}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb_path", type=Path, help="Prepared PDB to validate")
    parser.add_argument(
        "--custom-substructures",
        type=Path,
        default=None,
        help="Optional private OpenFF custom-substructures JSON proof of concept",
    )
    return parser.parse_args()


def main() -> None:
    """Run direct OpenFF validation."""
    args = parse_args()
    validate_openff(args.pdb_path, args.custom_substructures)


if __name__ == "__main__":
    main()
