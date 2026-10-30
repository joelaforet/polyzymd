"""Methods to cite for builds of modified proteins."""

from pathlib import Path

import yaml

from polyzymd.config.schema import SimulationConfig
from polyzymd.simulation.build_citations import cite

QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "config.yaml"
SAGE = {
    "force_field": ["openff-2.3.0.offxml"],
    "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
}


def _cited(protein: dict) -> str:
    data = yaml.safe_load(QUICKSTART.read_text())
    data["solutes"][0]["pdb_path"] = str(QUICKSTART.parent / "trpcage.pdb")
    data["solutes"][0].update(protein)
    return "\n".join(cite(SimulationConfig.model_validate(data), {}))


def test_an_unmodified_protein_cites_neither_mbuild_nor_pablo() -> None:
    text = _cited({})
    assert "mBuild" not in text and "openff-pablo" not in text
    assert "AshGC (OpenFF NAGL model" not in text


def test_an_octanoyl_site_cites_its_nagl_model_mbuild_and_pablo() -> None:
    text = _cited(
        {
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                    }
                }
            ],
            "modification_parameters": {"generic": SAGE},
        }
    )
    assert "AshGC (OpenFF NAGL model openff-gnn-am1bcc-1.0.0.pt)" in text
    assert "AM1-BCC model: I." in text
    assert "mBuild" in text and "openff-pablo" in text


def test_a_glycam_glycan_cites_glycam06(tmp_path) -> None:
    glycan = tmp_path / "glycam_G42666HT.pdb"
    glycan.write_text("REMARK placeholder\n")
    text = _cited(
        {
            "modifications": [
                {
                    "attach": {
                        "type": "glycan",
                        "site": "A:N1",
                        "atom_name": "ND2",
                        "reaction": "N-glycosylation",
                        "fragment": {"pdb_path": str(glycan)},
                    }
                }
            ],
            "modification_parameters": {
                "glycan": {"force_field": ["GLYCAM_06j-1.xml"], "charges": {"method": "library"}}
            },
        }
    )
    assert "GLYCAM06: a generalizable biomolecular force field" in text
    assert "PolyzyMD has no reference for this file" not in text


def test_a_mutation_alone_cites_mbuild_and_pablo() -> None:
    text = _cited({"mutations": [{"site": "A:K8", "to": "ALA"}]})
    assert "mBuild" in text and "openff-pablo" in text
    assert "AshGC (OpenFF NAGL model" not in text
