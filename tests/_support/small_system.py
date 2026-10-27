"""Build a small system that holds every built-in component type, for audit tests.

Trp-cage (the quickstart protein), a substrate SDF, polymer chains of two
monomers from pre-built charged SDFs, DMSO from PolyzyMD's bundled co-solvents,
water and NaCl. It builds in well under a minute, so tests can make it once and
then break copies of it on purpose.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"

#: Monomer A and B of the test polymer, written as 2-mers with explicit charges.
_POLYMER_SMILES = {
    "AA": "CC(C)(C(=O)OC)CC(C)(C)C(=O)OC",
    "AB": "CC(C)(C(=O)OC)CC(C)(C)C(=O)OCCO",
    "BB": "CC(C)(C(=O)OCCO)CC(C)(C)C(=O)OCCO",
}
SUBSTRATE_SMILES = "CC(=O)Oc1ccc(NC(C)=O)cc1"


def write_inputs(folder: Path) -> Path:
    """Write the protein, substrate, polymer SDFs and config into ``folder``; return the config."""
    from openff.toolkit import Molecule

    folder.mkdir(parents=True, exist_ok=True)
    (folder / "trpcage.pdb").write_text((QUICKSTART / "trpcage.pdb").read_text())

    substrate = Molecule.from_smiles(SUBSTRATE_SMILES)
    substrate.generate_conformers(n_conformers=1)
    # Place the substrate beside the protein, as a docked ligand would be.
    substrate.conformers[0] = substrate.conformers[0] + _protein_corner(folder / "trpcage.pdb")
    substrate.to_file(str(folder / "substrate.sdf"), file_format="sdf")

    polymers = folder / "polymers"
    polymers.mkdir(exist_ok=True)
    for sequence, smiles in _POLYMER_SMILES.items():
        chain = Molecule.from_smiles(smiles)
        chain.generate_conformers(n_conformers=1)
        chain.assign_partial_charges("gasteiger")
        chain.to_file(str(polymers / f"TP_seq={sequence}_2-mer_charged.sdf"), file_format="sdf")

    config: dict[str, Any] = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["name"] = "audit_small_system"
    sage = config["solvent"]["force_field"]
    config["solutes"] += [
        {
            "type": "small_molecule",
            "name": "substrate",
            "sdf_path": "substrate.sdf",
            "conformer_index": 0,
            "residue_name": "LIG",
            "force_field": sage,
            "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
        },
        {
            "type": "polymer",
            "type_prefix": "TP",
            "monomers": [
                {"label": "A", "probability": 0.5, "name": "MA"},
                {"label": "B", "probability": 0.5, "name": "HEMA"},
            ],
            "length": 2,
            "count": 6,
            "sdf_directory": "polymers",
            "packing": {"padding": 0.6, "tolerance": 2.0},
            "force_field": sage,
            "charges": {"method": "preset"},
        },
    ]
    config["solvent"]["box"]["padding"] = 0.8
    config["solvent"]["co_solvents"] = [
        {"name": "dmso", "concentration": 0.2, "charges": {"method": "preset"}}
    ]
    config["output"]["projects_directory"] = str(folder)
    config["output"]["scratch_directory"] = str(folder / "scratch")
    path = folder / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


def _protein_corner(pdb: Path):
    """Return an offset that puts a molecule just outside the protein's bounding box."""
    import numpy as np
    from openff.units import unit
    from openmm import app
    from openmm import unit as openmm_unit

    positions = app.PDBFile(str(pdb)).getPositions(asNumpy=True).value_in_unit(openmm_unit.angstrom)
    return (np.asarray(positions).max(axis=0) + 3.0) * unit.angstrom


def build(folder: Path, replicate: int = 1):
    """Build the small system; return the builder and its config."""
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.loader import load_config

    config = load_config(write_inputs(folder))
    builder = SystemBuilder.from_config(config)
    builder.build_from_config(config, working_dir=folder / "build", polymer_seed=replicate)
    return builder, config
