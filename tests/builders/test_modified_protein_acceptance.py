"""Acceptance builds of the mBuild protein demos and attach-movie systems.

Each case is a real config built with ``polyzymd build`` and checked by the audit.
They are slow and need files that are not in the repository (the
mbuild_protein_demos checkout and the attach-movie inputs), so they run only with
``POLYZYMD_ACCEPTANCE=1``. ``POLYZYMD_LAB`` names the folder that holds those
checkouts (default ``~/Shirts-Lab-Linux``); a case whose inputs are missing is
skipped. With ``POLYZYMD_ACCEPTANCE_RUN=1`` each build is also run briefly
(minimization and the config's short stages) through ``polyzymd run``, using the
Python named by ``POLYZYMD_SIM_PYTHON`` (an environment whose OpenMM matches the
GPU driver).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples" / "modified_proteins"
DATA = ROOT / "tests" / "data" / "modifying_proteins"
LAB = Path(os.environ.get("POLYZYMD_LAB", Path.home() / "Shirts-Lab-Linux"))
DEMOS = LAB / "mbuild_protein_demos"

SAGE = {
    "force_field": ["openff-2.3.0.offxml"],
    "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
}
GLYCAM = {"force_field": ["GLYCAM_06j-1.xml"], "charges": {"method": "library"}}
GLYCAN = {"type": "glycan", "fragment_atom_name": "C1", "fragment_resnum": 2}

#: The two FRET dyes of mbuild_protein_demos notebook 04.
DONOR = (
    "CC1(CCCCCC(=O)NCCC(=O)N2Cc3ccccc3C#Cc3ccccc32)C(/C=C/C=C2/N(CCCS(=O)(=O)[O-])"
    "c3ccc(S(=O)(=O)[O-])cc3C2(C)C)=[N+](CCCS(=O)(=O)[O-])c2ccc(S(=O)(=O)[O-])cc21"
)
ACCEPTOR = (
    "CC1(C)C(/C=C/C=C/C=C2\\N(CCCS(=O)(=O)[O-])c3ccc(S(=O)(=O)[O-])cc3C2(C)CCCCC(=O)"
    "NCCN2C(=O)C=CC2=O)=[N+](CCCS(=O)(=O)[O-])c2ccc(S(=O)(=O)[O-])cc21"
)

#: FN7-10 with the two mutations of notebook 03, which notebook 04 labels.
FN_MUTATIONS = [
    {"site": "A:S1381", "to": {"ccd_file": "structures/4II.cif"}},
    {"site": "A:S1500", "to": "CYS"},
]

#: The attach-movie LipA build: AzF and Cys mutations, the two FRET dyes, and a 22-unit
#: SBMA chain grafted on Lys44 through an NHS ester in its middle, with soluble chains.
LIPA = LAB / "mbuild_attach_movie" / "lipa" / "data" / "lipa_prepared.pdb"
NHS_AMINOLYSIS = "[N:1]([H:2])[H].[C:3](=[O:4])[O][N]>>[N:1]([H:2])[C:3]=[O:4]"
CAP = [{"label": "H", "smiles": "[*][H]", "probability": 1.0}]
SBMA_POLYMER = {
    "type": "polymer",
    "name": "pSBMA",
    "generation_mode": "dynamic",
    "type_prefix": "SBMA-EGPMA",
    "monomers": [
        {
            "label": "A",
            "name": "SBMA",
            "residue_name": "SBM",
            "smiles": "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]",
            "probability": 0.95,
        },
        {
            "label": "B",
            "name": "EGPMA",
            "residue_name": "EGP",
            "smiles": "[*]CC([*])(C)C(=O)OCCOc1ccccc1",
            "probability": 0.05,
        },
        {
            "label": "N",
            "name": "NHS",
            "residue_name": "NHS",
            "smiles": "[*]CC([*])(C)C(=O)ON1C(=O)CCC1=O",
            "probability": 0.0,
        },
    ],
    "end_groups": {"head": CAP, "tail": CAP},
    "length": 7,
    "count": 2,
    "force_field": ["openff-2.3.0.offxml"],
    "charges": SAGE["charges"],
}
#: Semaglutide: the C18 diacid linker (CCD KUT) on Lys20. Its CCD file comes
#: from the CCD cache openff-pablo keeps.
CCD_CACHE = Path(
    os.environ.get("POLYZYMD_CCD_CACHE", Path.home() / ".cache/openff-pablo/ccd_cache")
)
SEMAGLUTIDE = DEMOS / "semaglutide_apo.pdb"
BG505 = LAB / "Sprenger_Colab" / "wt-bg505"
BG505_TRIMER = (
    BG505
    / "glycosylated_mutants"
    / "structures_5FYL_minimized"
    / "trimer_with_fab_5FYL_assembled_ff14SB_minimized.pdb"
)
BG505_SITES = BG505 / "glycan_structures" / "gp140_glycans.csv"
BG505_GLYCANS = BG505 / "glycosylated_mutants" / "glycans_conect"
BG505_INPUTS = [BG505_TRIMER, BG505_SITES, BG505_GLYCANS]

CASES: dict[str, dict] = {
    # Notebook 03: canonical and non-canonical mutations, and a D form.
    "fn_mutations": {
        "pdb": DEMOS / "1fnf_prepared.pdb",
        "protein": {
            "mutations": [
                *FN_MUTATIONS,
                {"site": "A:S1367", "to": "LEU"},
                {"site": "A:K1275", "to": "ALA", "stereo": "D"},
            ],
            "modification_parameters": {"noncanonical": SAGE},
        },
    },
    # Notebook 04: the donor clicked onto AzF (merged), the acceptor on Cys.
    "fn_fret": {
        "pdb": DEMOS / "1fnf_prepared.pdb",
        "protein": {
            "mutations": FN_MUTATIONS,
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:4II1381",
                        "atom_name": "N3",
                        "reaction": "azide-alkyne triazole",
                        "merge": True,
                        "fragment": {"smiles": DONOR, "resname": "DON"},
                    }
                },
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:C1500",
                        "atom_name": "SG",
                        "reaction": "thiol-maleimide",
                        "fragment": {"smiles": ACCEPTOR, "resname": "ACC"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE, "noncanonical": SAGE},
        },
    },
    # The reaction movie: five named reactions on the mutated FN7-10.
    "fn_reactions": {
        "pdb": DEMOS / "1fnf_prepared.pdb",
        "glycans": [DEMOS / "glycans" / "glycam_G57321FI.pdb"],
        "protein": {
            "mutations": FN_MUTATIONS,
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:C1500",
                        "atom_name": "SG",
                        "reaction": "thiol-maleimide",
                        "fragment": {"smiles": "CCN1C(=O)C=CC1=O", "resname": "NEM"},
                    }
                },
                {"deprotonate": {"site": "A:K1275", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K1275",
                        "atom_name": "NZ",
                        "reaction": "amide coupling",
                        "fragment": {
                            "smiles": "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
                            "resname": "BTN",
                        },
                    }
                },
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:4II1381",
                        "atom_name": "N3",
                        "reaction": "azide-alkyne triazole",
                        "merge": True,
                        "fragment": {"smiles": "CC(=O)N1Cc2ccccc2C#Cc2ccccc21", "resname": "DBC"},
                    }
                },
                {
                    "attach": {
                        **GLYCAN,
                        "site": "A:N1183",
                        "atom_name": "ND2",
                        "reaction": "N-glycosylation",
                        "fragment": {"pdb_path": "structures/glycans/glycam_G57321FI.pdb"},
                    }
                },
                {
                    "attach": {
                        **GLYCAN,
                        "site": "A:S1367",
                        "atom_name": "OG",
                        "reaction": "O-glycosylation",
                        "fragment": {"pdb_path": "structures/glycans/glycam_G57321FI.pdb"},
                    }
                },
            ],
            "modification_parameters": {
                "generic": SAGE,
                "noncanonical": SAGE,
                "glycan": GLYCAM,
            },
        },
    },
    "lipa_fret_graft": {
        "pdb": LIPA,
        "solutes": [SBMA_POLYMER],
        "protein": {
            "mutations": [
                {"site": "A:N4", "to": {"ccd_file": "structures/4II.cif"}},
                {"site": "A:N174", "to": "CYS"},
            ],
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:4II4",
                        "atom_name": "N3",
                        "reaction": "azide-alkyne triazole",
                        "merge": True,
                        "fragment": {"smiles": DONOR, "resname": "DON"},
                    }
                },
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:C174",
                        "atom_name": "SG",
                        "reaction": "thiol-maleimide",
                        "fragment": {"smiles": ACCEPTOR, "resname": "ACC"},
                    }
                },
                {"deprotonate": {"site": "A:K44", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "polymer",
                        "site": "A:K44",
                        "atom_name": "NZ",
                        "reaction": NHS_AMINOLYSIS,
                        "attach_at": 11,
                        "fragment": {
                            "polymer": {
                                "monomers_from": "pSBMA",
                                "length": 22,
                                "sequence": {"fixed": {11: "NHS"}},
                            }
                        },
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE, "noncanonical": SAGE, "polymer": SAGE},
        },
    },
    # Semaglutide with Ala2 in place of Aib2 and the KUT linker on Lys20. mBuild's mutate
    # cannot place Aib: its CA carries two methyls, so there is no single side chain to swap.
    "semaglutide": {
        "pdb": SEMAGLUTIDE,
        "ccd": [CCD_CACHE / "KUT.cif"],
        "prepare": "semaglutide",
        "protein": {
            "modifications": [
                {"deprotonate": {"site": "A:K20", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K20",
                        "atom_name": "NZ",
                        "leaving_atom_names": ["HZ2"],
                        "fragment": {"ccd_file": "structures/KUT.cif", "link_atom": "C33"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE},
        },
    },
    # The BG505 SOSIP protomer (gp120 chain A, gp41 chain B) with its 20 glycans.
    "bg505_protomer": {"pdb": BG505_TRIMER, "glycan_inputs": BG505_INPUTS, "prepare": "bg505"},
}

#: The same LipA build relaxed after each change instead of once at the end.
CASES["lipa_fret_graft_stepwise"] = {
    **CASES["lipa_fret_graft"],
    "protein": {**CASES["lipa_fret_graft"]["protein"], "relax": "per_modification"},
}

#: The shipped examples, built as they are.
EXAMPLE_CASES = {
    "glycosylated_ubiquitin": "config_glycan.yaml",
    "acylated_ubiquitin": "config_octanoyl.yaml",
}

pytestmark = [
    pytest.mark.skipif(
        not os.environ.get("POLYZYMD_ACCEPTANCE"),
        reason="local acceptance build; set POLYZYMD_ACCEPTANCE=1",
    ),
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


def materialize(name: str, folder: Path) -> Path:
    """Write the config of case ``name`` (and its inputs) into ``folder``; return the config."""
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copytree(EXAMPLES / "structures", folder / "structures", dirs_exist_ok=True)
    if name in EXAMPLE_CASES:
        config = folder / EXAMPLE_CASES[name]
        shutil.copy(EXAMPLES / EXAMPLE_CASES[name], config)
        return config
    case = CASES[name]
    for path in [
        case["pdb"],
        *case.get("glycans", []),
        *case.get("glycan_inputs", []),
        *case.get("ccd", []),
    ]:
        if not Path(path).exists():
            pytest.skip(f"{path} is not on this machine (set POLYZYMD_LAB)")
    shutil.copy(DATA / "4II.cif", folder / "structures")
    for ccd in case.get("ccd", []):
        shutil.copy(ccd, folder / "structures")
    if case.get("prepare") == "bg505":
        case = {**case, "protein": _bg505_protomer(folder)}
    elif case.get("prepare") == "semaglutide":
        _ala2_semaglutide(case["pdb"], folder / "structures" / "protein.pdb")
    else:
        shutil.copy(case["pdb"], folder / "structures" / "protein.pdb")
    for glycan in case.get("glycans", []):
        shutil.copy(glycan, folder / "structures" / "glycans")
    config = yaml.safe_load((EXAMPLES / "config_glycan.yaml").read_text())
    config["name"] = name
    config["description"] = f"acceptance case {name}"
    protein = config["solutes"][0]
    for key in ("modifications", "modification_parameters"):
        protein.pop(key, None)
    protein.update({"name": name.replace("_", ""), "pdb_path": "structures/protein.pdb"})
    protein.update(case["protein"])
    config["solutes"].extend(case.get("solutes", []))
    path = folder / f"{name}.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


def _ala2_semaglutide(source: Path, target: Path) -> None:
    """Write semaglutide with Ala in place of Aib2 (mBuild cannot place Aib by mutation).

    Aib's CB1 methyl becomes Ala's CB; CB2 and its hydrogens go, and HA takes CB2's
    direction at 1.09 A from CA."""
    import numpy as np

    lines, ca, cb2 = [], None, None
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20] == "AIB":
            name = line[12:16].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if name == "CA":
                ca = xyz
            if name == "CB2":
                cb2 = xyz
            if name.startswith(("CB2", "HB2")):
                continue
            renamed = {"CB1": " CB ", "HB11": " HB1", "HB12": " HB2", "HB13": " HB3"}
            line = (
                "ATOM  "
                + line[6:12]
                + renamed.get(name, line[12:16])
                + line[16:17]
                + "ALA"
                + line[20:]
            )
        lines.append(line)
    ha = ca + 1.09 * (cb2 - ca) / np.linalg.norm(cb2 - ca)
    at = next(
        i
        for i, record in enumerate(lines)
        if record[17:20] == "ALA" and record[22:26].strip() == "2" and record[12:16].strip() == "CA"
    )
    template = lines[at]
    lines.insert(
        at + 1,
        f"{template[:12]} HA {template[16:30]}{ha[0]:8.3f}{ha[1]:8.3f}{ha[2]:8.3f}{template[54:76]} H",
    )
    serial = 0
    for i, record in enumerate(lines):
        if record.startswith(("ATOM", "HETATM")):
            serial += 1
            lines[i] = f"{record[:6]}{serial:5d}{record[11:]}"
    target.write_text("\n".join(lines) + "\n")


def _bg505_protomer(folder: Path) -> dict:
    """Write chains A and B of the BG505 trimer and their glycans; return the protein keys."""
    import csv

    keep, lines = set(), []
    for line in BG505_TRIMER.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[21] in "AB":
            keep.add(int(line[6:11]))
            lines.append(line)
        elif line.startswith("TER") and len(line) > 21 and line[21] in "AB":
            lines.append(line)
        elif line.startswith("CONECT"):
            numbers = [int(line[i : i + 5]) for i in range(6, len(line.rstrip()), 5)]
            if all(n in keep for n in numbers):
                lines.append(line)
    (folder / "structures" / "protein.pdb").write_text("\n".join(lines + ["END"]) + "\n")
    groups: dict[str, dict[str, list[str]]] = {}
    for row in csv.DictReader(open(BG505_SITES)):
        if int(row["residue_number_1"]) != 301:  # reverted to the wild-type T301
            groups.setdefault(row["segment_id_1"], {}).setdefault(
                row["glycan_accession"], []
            ).append(f"N{row['residue_number_1']}")
    for accession in {a for group in groups.values() for a in group}:
        shutil.copy(
            BG505_GLYCANS / f"glycam_{accession}.pdb",
            folder / "structures" / "glycans",
        )
    return {
        "modifications": [
            {
                "attach": {
                    **GLYCAN,
                    "reaction": "N-glycosylation",
                    "atom_name": "ND2",
                    "chains": [chain],
                    "glycans": glycans,
                }
            }
            for chain, glycans in sorted(groups.items())
        ],
        "modification_parameters": {"glycan": GLYCAM},
    }


def _cli(python: str, config: Path, *arguments: str) -> None:
    result = subprocess.run(
        [python, "-W", "ignore", "-m", "polyzymd.cli.main", *arguments, "-c", config.name],
        cwd=config.parent,
        env={**os.environ, "HOME": str(config.parent / "home")},
        capture_output=True,
        text=True,
        timeout=7200,
    )
    assert result.returncode == 0, f"{result.stdout[-4000:]}\n{result.stderr[-4000:]}"


@pytest.mark.parametrize("name", [*EXAMPLE_CASES, *CASES])
def test_the_case_builds_passes_the_audit_and_runs(tmp_path: Path, name: str) -> None:
    config = materialize(name, tmp_path / name)
    _cli(sys.executable, config, "build")
    (record_path,) = config.parent.glob("*_run1/build_provenance.json")
    audit = json.loads(record_path.read_text())["audit"]
    flagged = [f for f in audit["findings"] if f["status"] in ("FAIL", "WARN")]
    assert audit["passed"] and not flagged, flagged
    if os.environ.get("POLYZYMD_ACCEPTANCE_RUN"):
        _cli(os.environ.get("POLYZYMD_SIM_PYTHON", sys.executable), config, "run")
        (check,) = config.parent.glob("*_run1/production_0/production_0_system_check.json")
        assert json.loads(check.read_text())["matches_build"]


def test_the_bg505_trimer_plan_lists_63_glycans(tmp_path: Path) -> None:
    """The whole trimer as the glycosylation guide writes it: one entry for the three gp120
    chains, one for the three gp41 chains and one for the Fab light chains. ``validate
    --plans`` applies the 63 changes to the structure; nothing is built."""
    import csv

    for path in BG505_INPUTS:
        if not path.exists():
            pytest.skip(f"{path} is not on this machine (set POLYZYMD_LAB)")
    folder = tmp_path / "bg505_trimer"
    config = materialize("glycosylated_ubiquitin", folder)
    shutil.copy(BG505_TRIMER, folder / "structures" / "trimer.pdb")
    groups: dict[str, dict[str, list[str]]] = {"A": {}, "B": {}}
    for row in csv.DictReader(open(BG505_SITES)):
        if int(row["residue_number_1"]) != 301:  # reverted to the wild-type T301
            groups[row["segment_id_1"]].setdefault(row["glycan_accession"], []).append(
                f"N{row['residue_number_1']}"
            )
    for accession in {a for group in groups.values() for a in group} | {"G22768VO"}:
        shutil.copy(BG505_GLYCANS / f"glycam_{accession}.pdb", folder / "structures" / "glycans")
    link = {**GLYCAN, "reaction": "N-glycosylation", "atom_name": "ND2"}
    data = yaml.safe_load(config.read_text())
    protein = data["solutes"][0]
    protein["pdb_path"] = "structures/trimer.pdb"
    protein["modifications"] = [
        {"attach": {**link, "chains": ["A", "C", "E"], "glycans": groups["A"]}},
        {"attach": {**link, "chains": ["B", "D", "F"], "glycans": groups["B"]}},
        {"attach": {**link, "chains": ["I", "K", "M"], "glycans": {"G22768VO": ["N294"]}}},
    ]
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    result = subprocess.run(
        [
            sys.executable,
            "-W",
            "ignore",
            "-m",
            "polyzymd.cli.main",
            "validate",
            "-c",
            config.name,
            "--plans",
        ],
        cwd=folder,
        capture_output=True,
        text=True,
        timeout=7200,
    )
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]
    assert "(63, checked against its PDB)" in result.stdout
