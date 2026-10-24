"""The simulation environments (sim-cuda-*) have no MDAnalysis, RDKit or OpenFF.

A config with co-solvents must still load and pass the checks `run` and
`run-segment` make there.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import yaml

QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"


def test_the_vendored_selection_residues_match_mdanalysis() -> None:
    from MDAnalysis.core.selection import NucleicSelection, ProteinSelection

    from polyzymd.core.atom_groups import (
        NUCLEIC_SELECTION_RESIDUE_NAMES,
        PROTEIN_SELECTION_RESIDUE_NAMES,
    )

    assert PROTEIN_SELECTION_RESIDUE_NAMES == set(ProteinSelection.prot_res)
    assert NUCLEIC_SELECTION_RESIDUE_NAMES == set(NucleicSelection.nucl_res)


def test_a_cosolvent_config_loads_and_checks_without_mdanalysis(tmp_path) -> None:
    data = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    data["solvent"]["co_solvents"] = [
        {"name": "urea", "concentration": 2.0, "charges": {"method": "preset"}},
        {
            "name": "my_additive",
            "smiles": "CC(=O)NC",
            "count": 2,
            "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
        },
    ]
    (tmp_path / "trpcage.pdb").write_text((QUICKSTART / "trpcage.pdb").read_text())
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data))
    script = f"""
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name.split(".")[0] in ("MDAnalysis", "rdkit"):
            raise ModuleNotFoundError(f"No module named {{name!r}}")
sys.meta_path.insert(0, Block())
from polyzymd.config.loader import load_config
config = load_config({str(path)!r})
config.require_buildable("openmm", inputs=False)
print("ok", [c.residue_name for c in config.solvent.co_solvents])
"""
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0 and "ok" in result.stdout, result.stderr[-2000:]
