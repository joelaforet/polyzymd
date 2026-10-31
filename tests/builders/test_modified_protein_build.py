"""Proteins with mutations and modifications, from config to an audited build."""

from __future__ import annotations

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
QUICKSTART = ROOT / "examples" / "quickstart"
DATA = ROOT / "tests" / "data" / "modifying_proteins"
CHITOBIOSE = "glycam_G42666HT.pdb"

GLYCAM = {"force_field": ["GLYCAM_06j-1.xml"], "charges": {"method": "library"}}
SAGE = {
    "force_field": ["openff-2.3.0.offxml"],
    "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
}
N_GLYCAN = {
    "type": "glycan",
    "reaction": "N-glycosylation",
    "atom_name": "ND2",
    "fragment_atom_name": "C1",
    "fragment_resnum": 2,
}
O_GLYCAN = {**N_GLYCAN, "reaction": "O-glycosylation", "atom_name": "OG"}

#: Trp-cage Q5N with chitobiose on Asn5 (a glycans group) and on Ser13 (a site).
GLYCOSYLATED = {
    "mutations": [{"site": "A:Q5", "to": "ASN"}],
    "modifications": [
        {"attach": {**N_GLYCAN, "chains": ["A"], "glycans": {"G42666HT": ["N5"]}}},
        {
            "attach": {
                **O_GLYCAN,
                "sites": ["A:S13"],
                "fragment": {"pdb_path": f"structures/glycans/{CHITOBIOSE}"},
            }
        },
    ],
    "modification_parameters": {"glycan": GLYCAM},
}

pytestmark = [
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


def _folder(tmp_path: Path, protein: dict, sites_csv: str | None = None) -> Path:
    """A config folder: the quickstart with ``protein``'s keys added to the protein."""
    folder = tmp_path / "study"
    (folder / "structures" / "glycans").mkdir(parents=True)
    shutil.copy(QUICKSTART / "trpcage.pdb", folder)
    shutil.copy(DATA / CHITOBIOSE, folder / "structures" / "glycans")
    shutil.copy(DATA / "4II.cif", folder / "structures")
    config = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["solutes"][0].update(protein)
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    if sites_csv is not None:
        (folder / "sites.csv").write_text(sites_csv)
    return folder


def _build(folder: Path, *arguments: str) -> dict:
    """Run ``polyzymd build`` in ``folder`` and return its build record."""
    environment = {**os.environ, "HOME": str(folder / "home")}
    result = subprocess.run(
        [
            sys.executable,
            "-W",
            "ignore",
            "-m",
            "polyzymd.cli.main",
            "build",
            "-c",
            "config.yaml",
            *arguments,
        ],
        cwd=folder,
        env=environment,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"
    (record,) = folder.glob("*_run1/build_provenance.json")
    return json.loads(record.read_text())


def _passed(record: dict) -> None:
    audit = record["audit"]
    failed = [f for f in audit["findings"] if f["status"] in ("FAIL", "WARN")]
    assert audit["passed"] and not failed, failed


@pytest.mark.parametrize(
    "engine",
    [
        "openmm",
        pytest.param(
            "gromacs",
            marks=pytest.mark.skipif(
                shutil.which("gmx") is None, reason="GROMACS is not installed"
            ),
        ),
    ],
)
def test_a_glycosylated_protein_builds_and_passes_the_audit(tmp_path: Path, engine: str) -> None:
    folder = _folder(tmp_path, GLYCOSYLATED)
    record = _build(folder, *(("--format", "gromacs") if engine == "gromacs" else ()))
    _passed(record)
    modifications = record["provenance"]["protein_modifications"]
    assert [(d["kind"], d["site"]) for d in modifications["declared"]] == [
        ("mutate", "A:Q5"),
        ("attach", "A:N5"),
        ("attach", "A:S13"),
    ]
    assert [site["name"] for site in modifications["parameterization"]["glycam"]] == [
        "G001",
        "G002",
    ]
    checks = {f["check"]: f["status"] for f in record["audit"]["findings"]}
    assert {"parameters.one_four", "modifications.applied", "masses"} <= set(checks)
    # The build is compared with the config's own blocks, not just recorded.
    assert checks["modifications.declared"] == "PASS"
    (run,) = folder.glob("*_run1")
    if engine == "openmm":
        residues: dict[str, set[int]] = {}
        for line in (run / "solvated_system.pdb").read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                residues.setdefault(line[72:76].strip(), set()).add(int(line[22:26]))
        # Each attachment segment numbers its residues from 1; the protein keeps 1..20.
        assert residues["G001"] == residues["G002"] == {1, 2}
        assert residues["PROT"] == set(range(1, 21))
    else:
        protein_itp = max(run.glob("gromacs/*_MOL0.itp"))
        pairs, section = 0, None
        for line in protein_itp.read_text().splitlines():
            stripped = line.split(";")[0].strip()
            if stripped.startswith("["):
                section = stripped.strip("[] ").lower()
            elif section == "pairs" and stripped and stripped.split()[2] == "2":
                pairs += 1
        written = sum(site["n_pairs"] for site in modifications["parameterization"]["glycam"])
        assert pairs == written > 0
        # The .gro numbers G001's residues from 1, as the PDB does.
        groups, current = {}, None
        for line in (run / "gromacs" / "components.ndx").read_text().splitlines():
            if line.startswith("["):
                current = line.strip("[] ")
                groups[current] = []
            elif current:
                groups[current] += [int(n) for n in line.split()]
        gro = (max(run.glob("gromacs/*.gro"))).read_text().splitlines()[2:-1]
        g001 = {int(gro[n - 1][:5]) for n in groups["G001"]}
        assert g001 == {1, 2}


def test_an_acylated_lysine_builds_with_sage_and_nagl(tmp_path: Path) -> None:
    folder = _folder(
        tmp_path,
        {
            "modifications": [
                {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE},
        },
    )
    record = _build(folder)
    _passed(record)
    (site,) = record["provenance"]["protein_modifications"]["parameterization"]["sites"]
    assert site["name"] == "M001" and abs(site["shift_per_atom"]) < 0.005


#: NHS ester aminolysis: the amine bonds the ester carbon, N-hydroxysuccinimide leaves.
NHS_AMINOLYSIS = "[N:1]([H:2])[H].[C:3](=[O:4])[O][N]>>[N:1]([H:2])[C:3]=[O:4]"
CAP = [{"label": "H", "smiles": "[*][H]", "probability": 1.0}]
GRAFT_MONOMERS = [
    {
        "label": "A",
        "name": "SBMA",
        "residue_name": "SBM",
        "smiles": "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]",
        "probability": 1.0,
    },
    {
        "label": "N",
        "name": "NHS",
        "residue_name": "NHS",
        "smiles": "[*]CC([*])(C)C(=O)ON1C(=O)CCC1=O",
        "probability": 0.0,
    },
]
#: Trp-cage Lys8 carrying a 9-position SBMA chain through an NHS linker in the middle.
GRAFTED = {
    "modifications": [
        {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
        {
            "attach": {
                "type": "polymer",
                "site": "A:K8",
                "atom_name": "NZ",
                "reaction": NHS_AMINOLYSIS,
                "attach_at": 5,
                "fragment": {
                    "polymer": {
                        "monomers": GRAFT_MONOMERS,
                        "end_groups": {"head": CAP, "tail": CAP},
                        "length": 9,
                        "sequence": {"fixed": {5: "NHS"}},
                        "random_seed": 3,
                    }
                },
            }
        },
    ],
    "modification_parameters": {"polymer": SAGE},
}


def test_a_grafted_polymer_chain_builds_as_one_site_and_passes_the_audit(tmp_path: Path) -> None:
    record = _build(_folder(tmp_path, GRAFTED))
    _passed(record)
    modifications = record["provenance"]["protein_modifications"]
    (site,) = modifications["parameterization"]["sites"]
    assert site["name"] == "P001" and abs(site["shift_per_atom"]) < 0.005
    (graft,) = [e["attach"] for e in modifications["applied"]["modifications"] if "attach" in e]
    assert graft["fragment"]["polymer"]["drawn"] == "HAAANAAAH"
    assert graft["attach_at"] == 5
    declared = {f["check"]: f for f in record["audit"]["findings"]}["modifications.declared"]
    assert declared["status"] == "PASS", declared


def test_a_click_onto_a_noncanonical_residue_merges_into_one_site(tmp_path: Path) -> None:
    folder = _folder(
        tmp_path,
        {
            "mutations": [{"site": "A:Y3", "to": {"ccd_file": "structures/4II.cif"}}],
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:4II3",
                        "atom_name": "N3",
                        "reaction": "azide-alkyne triazole",
                        "merge": True,
                        "fragment": {"smiles": "C#CCO", "resname": "PRP"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE, "noncanonical": SAGE},
        },
    )
    record = _build(folder)
    _passed(record)
    (site,) = record["provenance"]["protein_modifications"]["parameterization"]["sites"]
    assert site["type"] == "noncanonical"


def test_a_sites_file_declares_the_same_modifications_as_the_inline_form(tmp_path: Path) -> None:
    from polyzymd.config.loader import load_config
    from polyzymd.config.modifications import declared_modifications
    from polyzymd.modifying_proteins import expand

    inline = load_config(_folder(tmp_path / "inline", GLYCOSYLATED) / "config.yaml").enzyme
    from_file = {
        **GLYCOSYLATED,
        "modifications": [
            {"attach": {**N_GLYCAN, "sites_file": "sites.csv"}},
            GLYCOSYLATED["modifications"][1],
        ],
    }
    folder = _folder(tmp_path / "csv", from_file, sites_csv="site,glycan\nA:N5,G42666HT\n")
    tabled = load_config(folder / "config.yaml").enzyme

    def declared(solute) -> list[tuple]:
        return [
            (d["kind"], d["site"], (d.get("fragment") or {}).get("sha256"))
            for d in declared_modifications(
                solute.mutations, solute.modifications, solute.modification_parameters
            )
        ]

    from polyzymd.builders.modified_protein import config_spec

    def expanded(solute) -> list[tuple]:
        _, modifications = expand(config_spec(solute))
        return [
            (verb, args["site"], Path(args["fragment"]["pdb_path"]).name)
            for entry in modifications
            for verb, args in entry.items()
        ]

    assert declared(inline) == declared(tabled)
    assert expanded(inline) == expanded(tabled)


@pytest.mark.skipif(
    not os.environ.get("POLYZYMD_ACCEPTANCE"),
    reason="local acceptance build; set POLYZYMD_ACCEPTANCE=1 (needs the demos and BG505 files)",
)
def test_ubiquitin_with_a_high_mannose_glycan(tmp_path: Path) -> None:
    lab = Path.home() / "Shirts-Lab-Linux"
    ubiquitin = lab / "mbuild_protein_demos" / "1ubq_protonated.pdb"
    man5 = lab / "Sprenger_Colab/wt-bg505/glycosylated_mutants/glycans_conect/glycam_G55220VL.pdb"
    if not (ubiquitin.is_file() and man5.is_file()):
        pytest.skip("acceptance inputs are not on this machine")
    folder = _folder(
        tmp_path,
        {
            "pdb_path": "1ubq_protonated.pdb",
            "modifications": [
                {"attach": {**N_GLYCAN, "glycans": {"G55220VL": ["A:N60"]}}},
            ],
            "modification_parameters": {"glycan": GLYCAM},
        },
    )
    shutil.copy(ubiquitin, folder)
    shutil.copy(man5, folder / "structures" / "glycans")
    record = _build(folder)
    _passed(record)


@pytest.mark.skipif(
    not os.environ.get("POLYZYMD_ACCEPTANCE"),
    reason="local acceptance build; set POLYZYMD_ACCEPTANCE=1 (needs the LipA movie files)",
)
def test_lipa_with_a_grafted_20mer_and_soluble_chains_of_the_same_polymer(tmp_path: Path) -> None:
    """The LipA graft of the attach movie: a 20-mer on Lys44 through an NHS linker in its middle.

    The chain is drawn for replicate 1 (no random_seed)."""
    lipa = Path.home() / "Shirts-Lab-Linux/mbuild_attach_movie/lipa/data/lipa_prepared.pdb"
    if not lipa.is_file():
        pytest.skip("acceptance inputs are not on this machine")
    egpma = {
        "label": "B",
        "name": "EGPMA",
        "residue_name": "EGP",
        "smiles": "[*]CC([*])(C)C(=O)OCCOc1ccccc1",
        "probability": 0.05,
    }
    monomers = [{**GRAFT_MONOMERS[0], "probability": 0.95}, egpma, GRAFT_MONOMERS[1]]
    graft = copy.deepcopy(GRAFTED)
    attach = graft["modifications"][1]["attach"]
    attach.update(site="A:K44", attach_at=11)
    attach["fragment"] = {
        "polymer": {
            "monomers_from": "pSBMA",
            "length": 22,
            "sequence": {"fixed": {11: "NHS"}},
        }
    }
    graft["modifications"][0]["deprotonate"]["site"] = "A:K44"
    folder = _folder(tmp_path, {"pdb_path": "lipa_prepared.pdb", **graft})
    shutil.copy(lipa, folder)
    config = yaml.safe_load((folder / "config.yaml").read_text())
    config["solutes"].append(
        {
            "type": "polymer",
            "name": "pSBMA",
            "generation_mode": "dynamic",
            "type_prefix": "SBMA-EGPMA",
            "monomers": monomers,
            "end_groups": {"head": CAP, "tail": CAP},
            "length": 7,
            "count": 2,
            "force_field": ["openff-2.3.0.offxml"],
            "charges": SAGE["charges"],
        }
    )
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    record = _build(folder)
    _passed(record)
    (site,) = record["provenance"]["protein_modifications"]["parameterization"]["sites"]
    assert site["name"] == "P001"
    findings = {f["check"]: f["detail"] for f in record["audit"]["findings"]}
    assert "grafted in 1 chain" in findings["composition.polymer_mass_fraction"]
    assert "grafted P001" in findings["composition.polymer_sequences"]


def _validate(folder: Path, *arguments: str) -> subprocess.CompletedProcess:
    environment = {**os.environ, "HOME": str(folder / "home")}
    return subprocess.run(
        [
            sys.executable,
            "-W",
            "ignore",
            "-m",
            "polyzymd.cli.main",
            "validate",
            "-c",
            "config.yaml",
            "--plans",
            *arguments,
        ],
        cwd=folder,
        env=environment,
        capture_output=True,
        text=True,
        timeout=600,
    )


def test_validate_plans_lists_one_row_per_site(tmp_path: Path) -> None:
    folder = _folder(tmp_path, GLYCOSYLATED)
    text = _validate(folder)
    assert text.returncode == 0, text.stderr
    assert "Modifications of trpcage (3, checked against its PDB)" in text.stdout
    assert "G002  attach      A:S13" in text.stdout
    agent = _validate(folder, "--format", "agent")
    rows = json.loads(agent.stdout)["modifications"]["trpcage"]
    assert [(row["segment"], row["site"], row["atoms_removed"]) for row in rows[1:]] == [
        ("G001", "A:N5", ["HD22"]),
        ("G002", "A:S13", ["HG"]),
    ]
    assert rows[1]["parameters"]["force_field"] == ["GLYCAM_06j-1.xml"]


def test_validate_plans_names_a_site_that_is_not_in_the_structure(tmp_path: Path) -> None:
    wrong = copy.deepcopy(GLYCOSYLATED)
    wrong["modifications"][0]["attach"]["glycans"] = {"G42666HT": ["N6"]}
    result = _validate(_folder(tmp_path, wrong))
    assert result.returncode != 0
    assert "residue 6 of chain A is TRP" in result.stdout + result.stderr


def test_a_graft_without_random_seed_is_seeded_by_the_replicate() -> None:
    from polyzymd.builders.modified_protein import config_spec
    from polyzymd.config.schema import ProteinSolute

    def chain(protein: dict, replicate: int) -> dict:
        solute = ProteinSolute.model_validate(
            {
                "type": "protein",
                "pdb_path": str(QUICKSTART / "trpcage.pdb"),
                "force_field": ["ff14sb_off_impropers_0.0.4.offxml"],
                "charges": {"method": "library"},
                **protein,
            }
        )
        (graft,) = [
            e["attach"] for e in config_spec(solute, replicate)["modifications"] if "attach" in e
        ]
        return graft["fragment"]["polymer"]

    per_replicate = copy.deepcopy(GRAFTED)
    del per_replicate["modifications"][1]["attach"]["fragment"]["polymer"]["random_seed"]
    assert chain(per_replicate, 2)["replicate"] == 2
    assert "replicate" not in chain(GRAFTED, 2)


#: The donor dye of the FRET labeling demo (mbuild_protein_demos notebook 04): a cyanine
#: with a DBCO alkyne for the strain-promoted click onto p-azido-phenylalanine (4II).
DONOR = (
    "CC1(CCCCCC(=O)NCCC(=O)N2Cc3ccccc3C#Cc3ccccc32)C(/C=C/C=C2/N(CCCS(=O)(=O)[O-])"
    "c3ccc(S(=O)(=O)[O-])cc3C2(C)C)=[N+](CCCS(=O)(=O)[O-])c2ccc(S(=O)(=O)[O-])cc21"
)


def test_the_fret_donor_clicks_onto_azidophenylalanine(tmp_path: Path) -> None:
    folder = _folder(
        tmp_path,
        {
            "mutations": [{"site": "A:Y3", "to": {"ccd_file": "structures/4II.cif"}}],
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:4II3",
                        "atom_name": "N3",
                        "reaction": "azide-alkyne triazole",
                        "merge": True,
                        "fragment": {"smiles": DONOR, "resname": "DON"},
                    }
                }
            ],
            "modification_parameters": {"generic": SAGE, "noncanonical": SAGE},
        },
    )
    record = _build(folder)
    _passed(record)
    (site,) = record["provenance"]["protein_modifications"]["parameterization"]["sites"]
    assert site["type"] == "noncanonical" and abs(site["shift_per_atom"]) < 0.005


def _two_chains(folder: Path) -> None:
    """Replace the folder's trpcage.pdb with two copies: chain A, and chain B 3 nm away."""
    atoms = [
        line
        for line in (QUICKSTART / "trpcage.pdb").read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    shifted = [
        f"{line[:21]}B{line[22:30]}{float(line[30:38]) + 30:8.3f}{line[38:]}" for line in atoms
    ]
    (folder / "trpcage.pdb").write_text(
        "\n".join(atoms + ["TER"] + shifted + ["TER", "END"]) + "\n"
    )


@pytest.mark.parametrize(
    "engine",
    [
        "openmm",
        pytest.param(
            "gromacs",
            marks=pytest.mark.skipif(
                shutil.which("gmx") is None, reason="GROMACS is not installed"
            ),
        ),
    ],
)
def test_a_two_chain_protein_builds_one_molecule_per_chain(tmp_path: Path, engine: str) -> None:
    folder = _folder(
        tmp_path,
        {
            "mutations": [{"site": "B:Q5", "to": "ASN"}],
            "modifications": [
                {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                    }
                },
                {"attach": {**N_GLYCAN, "chains": ["B"], "glycans": {"G42666HT": ["N5"]}}},
            ],
            "modification_parameters": {"generic": SAGE, "glycan": GLYCAM},
        },
    )
    _two_chains(folder)
    record = _build(folder, *(("--format", "gromacs") if engine == "gromacs" else ()))
    _passed(record)
    modifications = record["provenance"]["protein_modifications"]
    first, second = modifications["parameterization"]["molecules"]
    assert first[0] == 0 and first[1] == second[0]
    (protein,) = [c for c in record["components"] if c["type"] == "protein"]
    assert protein["n_molecules"] == 2


def test_validate_plans_reports_a_reaction_that_does_not_apply(tmp_path: Path) -> None:
    folder = _folder(
        tmp_path,
        {
            "modifications": [
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:S13",
                        "atom_name": "OG",
                        "reaction": "thiol-maleimide",
                        "fragment": {"smiles": "CCN1C(=O)C=CC1=O", "resname": "NEM"},
                    }
                }
            ],
            "modification_parameters": {"generic": SAGE},
        },
    )
    result = _validate(folder)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "error:" in output and "fix:" in output and "Traceback" not in output


def test_validate_format_agent_needs_plans(tmp_path: Path) -> None:
    folder = _folder(tmp_path, {})
    environment = {**os.environ, "HOME": str(folder / "home")}
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "polyzymd.cli.main",
            "validate",
            "-c",
            "config.yaml",
            "--format",
            "agent",
        ],
        cwd=folder,
        env=environment,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode != 0 and "use it with --plans" in result.stderr


def test_the_audit_fails_a_build_that_drops_a_declared_glycan(tmp_path, monkeypatch) -> None:
    from click.testing import CliRunner

    import polyzymd.builders.modified_protein as modified_protein
    from polyzymd.cli.main import cli

    folder = _folder(tmp_path, GLYCOSYLATED)
    original = modified_protein.config_spec

    def without_the_last_glycan(solute, *args, **kwargs):
        spec = original(solute, *args, **kwargs)
        spec["modifications"] = spec["modifications"][:-1]
        return spec

    monkeypatch.setattr(modified_protein, "config_spec", without_the_last_glycan)
    monkeypatch.chdir(folder)
    monkeypatch.setenv("HOME", str(folder / "home"))
    result = CliRunner().invoke(cli, ["build", "-c", "config.yaml"])
    assert result.exit_code != 0
    (failed,) = folder.glob("*_run1/build_provenance.failed.json")
    findings = json.loads(failed.read_text())["audit"]["findings"]
    declared = [f for f in findings if f["check"] == "modifications.declared"]
    assert declared and declared[0]["status"] == "FAIL", declared


@pytest.mark.skipif(
    not os.environ.get("POLYZYMD_ACCEPTANCE"),
    reason="local acceptance build; set POLYZYMD_ACCEPTANCE=1 (needs the BG505 files)",
)
def test_a_bg505_protomer_with_its_twenty_glycans(tmp_path: Path) -> None:
    """gp120 + gp41 of BG505 SOSIP (chains A and B, one molecule through the SOS
    disulfide) with the 20 glycans of the Sprenger build, on GLYCAM06j."""
    import csv

    bg505 = Path.home() / "Shirts-Lab-Linux" / "Sprenger_Colab" / "wt-bg505"
    trimer = (
        bg505
        / "glycosylated_mutants"
        / "structures_5FYL_minimized"
        / "trimer_with_fab_5FYL_assembled_ff14SB_minimized.pdb"
    )
    if not trimer.is_file():
        pytest.skip("BG505 inputs are not on this machine")
    folder = _folder(tmp_path, {})
    keep, lines = set(), []
    for line in trimer.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[21] in "AB":
            keep.add(int(line[6:11]))
            lines.append(line)
        elif line.startswith("TER") and len(line) > 21 and line[21] in "AB":
            lines.append(line)
        elif line.startswith("CONECT"):
            numbers = [int(line[i : i + 5]) for i in range(6, len(line.rstrip()), 5)]
            if all(n in keep for n in numbers):
                lines.append(line)
    (folder / "protomer.pdb").write_text("\n".join(lines + ["END"]) + "\n")
    groups: dict[str, dict[str, list[str]]] = {}
    for row in csv.DictReader(open(bg505 / "glycan_structures" / "gp140_glycans.csv")):
        if int(row["residue_number_1"]) != 301:  # reverted to the wild-type T301
            groups.setdefault(row["segment_id_1"], {}).setdefault(
                row["glycan_accession"], []
            ).append(f"N{row['residue_number_1']}")
    library = folder / "structures" / "glycans"
    for accession in {a for group in groups.values() for a in group}:
        shutil.copy(
            bg505 / "glycosylated_mutants" / "glycans_conect" / f"glycam_{accession}.pdb", library
        )
    config = yaml.safe_load((folder / "config.yaml").read_text())
    config["solutes"][0].update(
        {
            "pdb_path": "protomer.pdb",
            "modifications": [
                {"attach": {**N_GLYCAN, "chains": [chain], "glycans": glycans}}
                for chain, glycans in sorted(groups.items())
            ],
            "modification_parameters": {"glycan": GLYCAM},
        }
    )
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    record = _build(folder)
    _passed(record)
    assert len(record["provenance"]["protein_modifications"]["parameterization"]["glycam"]) == 20
    # Each glycan is placed clear of the protein and the glycans placed before it.
    statuses = {f["check"]: f["status"] for f in record["audit"]["findings"]}
    assert statuses["geometry.modified_protein_clash"] == "PASS"


def test_adjacent_one_residue_glycans_stay_separate_residues(tmp_path: Path) -> None:
    """Three one-residue glycans in a row are each "resid 1" of their segment; the PDB
    separates them with TER so a reader does not merge them as alternate locations."""
    glcnac = {"pdb_path": "structures/glycans/glycam_G49108TO.pdb"}
    folder = _folder(
        tmp_path,
        {
            "mutations": [{"site": "A:Q5", "to": "ASN"}],
            "modifications": [
                {"attach": {**N_GLYCAN, "site": "A:N5", "fragment": glcnac}},
                {"attach": {**O_GLYCAN, "sites": ["A:S13", "A:S14"], "fragment": glcnac}},
            ],
            "modification_parameters": {"glycan": GLYCAM},
        },
    )
    shutil.copy(DATA / "glycam_G49108TO.pdb", folder / "structures" / "glycans")
    record = _build(folder)
    _passed(record)
    (run,) = folder.glob("*_run1")
    import openmm.app as app

    residues = [
        (r.chain.index, r.id, r.name)
        for r in app.PDBFile(str(run / "solvated_system.pdb")).topology.residues()
        if r.name == "0YB"
    ]
    assert len(residues) == 3 and len(set(residues)) == 3


def test_the_audit_fails_a_build_whose_expansion_drops_a_group_site(tmp_path, monkeypatch) -> None:
    """The audit's declared side is the config's own expansion, so a wrong expansion
    in the builder (here: the glycans group loses its site) fails modifications.declared."""
    from click.testing import CliRunner

    import polyzymd.modifying_proteins.spec as spec
    from polyzymd.cli.main import cli

    folder = _folder(tmp_path, GLYCOSYLATED)
    original = spec.expand

    def losing_a_group_site(blocks, base_dir="."):
        mutations, modifications = original(blocks, base_dir)
        grouped = any(
            "chains" in (entry.get("attach") or {}) for entry in blocks.get("modifications") or []
        )
        return mutations, (modifications[1:] if grouped else modifications)

    monkeypatch.setattr(spec, "expand", losing_a_group_site)
    monkeypatch.chdir(folder)
    monkeypatch.setenv("HOME", str(folder / "home"))
    result = CliRunner().invoke(cli, ["build", "-c", "config.yaml"])
    assert result.exit_code != 0
    (failed,) = folder.glob("*_run1/build_provenance.failed.json")
    findings = json.loads(failed.read_text())["audit"]["findings"]
    declared = [f for f in findings if f["check"] == "modifications.declared"]
    assert declared and declared[0]["status"] == "FAIL", declared


def _donor_click(relax: str) -> dict:
    return {
        "relax": relax,
        "mutations": [{"site": "A:Y3", "to": {"ccd_file": "structures/4II.cif"}}],
        "modifications": [
            {
                "attach": {
                    "type": "generic",
                    "site": "A:4II3",
                    "atom_name": "N3",
                    "reaction": "azide-alkyne triazole",
                    "merge": True,
                    "fragment": {"smiles": DONOR, "resname": "DON"},
                }
            }
        ],
        "modification_parameters": {"generic": SAGE, "noncanonical": SAGE},
    }


def test_joint_relaxation_places_the_dye_clear_of_the_protein(tmp_path: Path) -> None:
    record = _build(_folder(tmp_path, _donor_click("joint")))
    _passed(record)
    statuses = {f["check"]: f["status"] for f in record["audit"]["findings"]}
    assert statuses["geometry.modified_protein_clash"] == "PASS"


def test_the_audit_fails_a_build_with_a_hydrogen_placed_through_the_protein(
    tmp_path, monkeypatch
) -> None:
    """A fragment hydrogen 0.6 Å from a protein hydrogen, every heavy atom clear,
    fails the build's clash check."""
    import numpy as np
    from click.testing import CliRunner

    import polyzymd.builders.modified_protein as modified_protein
    from polyzymd.cli.main import cli

    original = modified_protein.apply_modifications

    def with_a_hydrogen_overlap(solute, base_dir, relax, *args, **kwargs):
        protein = original(solute, base_dir, relax, *args, **kwargs)
        if relax == "none":
            return protein  # the audit's replay
        atoms = list(protein.mbuild.particles())
        hydrogen = next(a for a in atoms if a.parent.name == "OC8" and a.name.startswith("H"))
        target = next(a for a in atoms if a.parent.name == "TRP" and a.name == "HE1")
        (parent,) = target.direct_bonds()
        away = (target.pos - parent.pos) / np.linalg.norm(target.pos - parent.pos)
        hydrogen.pos = target.pos + 0.06 * away
        return protein

    monkeypatch.setattr(modified_protein, "apply_modifications", with_a_hydrogen_overlap)
    folder = _folder(
        tmp_path,
        {
            "modifications": [
                {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE},
        },
    )
    monkeypatch.chdir(folder)
    monkeypatch.setenv("HOME", str(folder / "home"))
    result = CliRunner().invoke(cli, ["build", "-c", "config.yaml"])
    assert result.exit_code != 0
    (failed,) = folder.glob("*_run1/build_provenance.failed.json")
    findings = json.loads(failed.read_text())["audit"]["findings"]
    (clash,) = [f for f in findings if f["check"] == "geometry.modified_protein_clash"]
    assert clash["status"] == "FAIL", clash
    assert "1 hydrogen-hydrogen pairs closer than 0.8 Å" in clash["detail"]
    assert "heavy-heavy" not in clash["detail"] and "TRP6:HE1" in clash["detail"]


def test_joint_relaxation_places_each_change_with_mbuilds_torsion_search(
    tmp_path, monkeypatch
) -> None:
    """Joint mode keeps mBuild's clash-avoiding placement (relax) and skips only the
    per-change minimization, then relaxes everything once."""
    from polyzymd.builders.modified_protein import apply_modifications
    from polyzymd.config.loader import load_config
    from polyzymd.modifying_proteins import api

    calls = []
    for name in ("mutate", "attach"):
        original = getattr(api, name)

        def spy(*args, _original=original, _name=name, **kwargs):
            calls.append((_name, kwargs.get("relax"), kwargs.get("minimize")))
            return _original(*args, **kwargs)

        monkeypatch.setattr(api, name, spy)
    monkeypatch.setattr(api, "relax_all", lambda protein: calls.append(("relax_all", None, None)))
    folder = _folder(tmp_path, GLYCOSYLATED)
    solute = load_config(folder / "config.yaml").enzyme
    apply_modifications(solute, folder, "joint")
    assert calls == [
        ("mutate", True, False),
        ("attach", True, False),
        ("attach", True, False),
        ("relax_all", None, None),
    ]


def test_sites_use_input_numbers_and_restraints_the_written_ones(tmp_path: Path) -> None:
    """A protein whose PDB starts at residue 2: modification sites name the input's
    residue numbers (A:K9), while restraints resolve against the written structure,
    which numbers the protein from 1 as PolyzyMD 1.3 did (so input residue 2 is resid 1)."""
    folder = _folder(
        tmp_path,
        {
            "modifications": [
                {"deprotonate": {"site": "A:K9", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K9",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)C", "resname": "ACE"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE},
        },
    )
    lines = (QUICKSTART / "trpcage.pdb").read_text().splitlines()
    shifted = [
        (
            f"{line[:22]}{int(line[22:26]) + 1:4d}{line[26:]}"
            if line.startswith(("ATOM", "HETATM"))
            else line
        )
        for line in lines
    ]
    (folder / "trpcage.pdb").write_text("\n".join(shifted) + "\n")
    config = yaml.safe_load((folder / "config.yaml").read_text())
    config["restraints"] = [
        {
            "type": "flat_bottom",
            "name": "ends",
            "atom1": {"component": "trpcage", "selection": "protein and resid 1 and name CA"},
            "atom2": {"component": "trpcage", "selection": "protein and resid 20 and name CA"},
            "distance": 30.0,
        }
    ]
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    record = _build(folder)
    _passed(record)
    statuses = {f["check"]: f["status"] for f in record["audit"]["findings"]}
    assert statuses["restraint.atom_identity"] == "PASS"
    (run,) = folder.glob("*_run1")
    written = {
        int(line[22:26]): line[17:20]
        for line in (run / "solvated_system.pdb").read_text().splitlines()
        if line.startswith("ATOM") and line[72:76].strip() == "PROT"
    }
    modified = {
        int(line[22:26]): line[17:20]
        for line in (run / "modified_protein.pdb").read_text().splitlines()
        if line.startswith("ATOM")
    }
    assert written[1] == modified[2] == "ASN" and written[8] == modified[9] == "LYS"


def test_a_fragment_from_a_local_ccd_file_builds(tmp_path: Path) -> None:
    """{ccd_file, link_atom} fragments (a component mBuild does not bundle, such as a
    linker) are read from the file and recorded with its hash."""
    folder = _folder(
        tmp_path,
        {
            "modifications": [
                {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"ccd_file": "structures/4II.cif", "link_atom": "CB"},
                    }
                },
            ],
            "modification_parameters": {"generic": SAGE},
        },
    )
    record = _build(folder)
    _passed(record)
    statuses = {f["check"]: f["status"] for f in record["audit"]["findings"]}
    assert statuses["modifications.declared"] == "PASS"
    attach = record["provenance"]["protein_modifications"]["declared"][1]
    assert attach["fragment"]["link_atom"] == "CB" and len(attach["fragment"]["sha256"]) == 64


def test_the_audit_rebuild_does_not_log_placement_warnings(tmp_path: Path) -> None:
    """The audit replays the modifications unrelaxed; mBuild's warnings about close
    contacts in that replay do not describe the built structure and stay out of the log."""
    folder = _folder(tmp_path, GLYCOSYLATED)
    environment = {**os.environ, "HOME": str(folder / "home")}
    result = subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "polyzymd.cli.main", "build", "-c", "config.yaml"],
        cwd=folder,
        env=environment,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    assert result.returncode == 0, result.stderr
    log = result.stdout + result.stderr
    marker = "The audit rebuilds the modified protein without relaxation"
    assert marker in log
    after = log[log.index(marker) :]
    assert "within 2.0 A" not in after and "Relax the structure" not in after
