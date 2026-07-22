"""Opt-in final mixed glycan/polymer conjugation acceptance test."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml

from tests._support.conjugation_fast_e2e import FIXTURE_CONFIG, FIXTURE_DIR

GROMACS_E2E_TIMEOUT_SECONDS = 300
GLYCAM_RESIDUE_NAMES = frozenset({"0YB", "4YB", "ROH", "NLN", "NAG", "BMA", "MAN"})
POLYMER_RESIDUE_NAMES = frozenset({"LYX", "NHX", "SBM", "SB1", "EGP", "EG1", "NHS"})
COSOLVENT_RESIDUE_NAMES = frozenset({"EOH", "ETH", "ETOH", "UNL"})
WATER_RESIDUE_NAMES = frozenset({"HOH", "WAT", "SOL"})
ION_RESIDUE_NAMES = frozenset({"NA", "CL", "Na+", "Cl-"})


@pytest.mark.slow
@pytest.mark.final_e2e
def test_final_1ubq_mixed_glycan_polymer_e2e(tmp_path: Path) -> None:
    """Run the final opt-in 1UBQ mixed GLYCAM/OpenFF acceptance workflow."""

    if os.environ.get("POLYZYMD_RUN_FINAL_E2E") != "1":
        pytest.skip("Set POLYZYMD_RUN_FINAL_E2E=1 to run the final E2E acceptance test")
    _require_final_e2e_prerequisites()

    from polyzymd.builders.conjugation.api import build_conjugate_from_config
    from polyzymd.builders.conjugation.pablo.parameterization import (
        InterchangeParameterizationSettings,
    )
    from polyzymd.builders.conjugation.relaxation import ConjugateRelaxationSettings
    from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings
    from polyzymd.config.loader import load_config
    from polyzymd.exporters.gromacs import GromacsExporter

    config_path = _write_final_runtime_config(tmp_path)
    config = load_config(config_path)
    settings = ConjugatedPolymerSystemSettings(
        create_final_interchange=True,
        force_regenerate_conjugate_polymer=True,
        conjugate_polymerist_max_retries=3,
        conjugate_polymerist_energy_minimize=True,
        conjugate_parameterization=InterchangeParameterizationSettings(
            small_molecule_force_field="openff-2.3.0.offxml",
        ),
        relaxation=ConjugateRelaxationSettings(platform_name="CPU"),
    )

    result = build_conjugate_from_config(
        config,
        output_dir=tmp_path / "build",
        settings=settings,
        free_polymer_seed=202610,
    )
    bundle = result.exact_export_bundle
    if bundle is None:
        raise AssertionError("Final E2E expected an exact mixed-overlay OpenMM bundle")

    openmm_summary = _run_short_openmm_exact_bundle(bundle, tmp_path / "openmm")
    exported = GromacsExporter(bundle, config).export(
        tmp_path / "gromacs",
        prefix="final_e2e_1ubq",
        gmx_command=shutil.which("gmx") or "gmx",
    )
    gromacs_summary = _run_gromacs_acceptance(exported, tmp_path / "gromacs")
    summary = _final_e2e_summary(result, exported, openmm_summary, gromacs_summary)
    summary_path = tmp_path / "final_e2e_summary.json"
    summary["artifacts"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    _assert_final_e2e_summary(summary)
    assert summary["components"]["attached_glycans"] == 2
    assert summary["components"]["attached_polymers"] == 2
    assert summary["components"]["free_polymers"] == 1
    assert summary["components"]["cosolvents"] >= 1
    assert summary["components"]["water_atoms"] > 0
    assert summary["components"]["ion_atoms"] > 0
    assert summary["ownership"]["glycam_atoms"] > 0
    assert summary["parity"]["conflict_count"] == 0
    assert summary["exact_gromacs_audit"]["pair_mismatch_count"] == 0
    assert summary["exact_gromacs_audit"]["exclusion_mismatch_count"] == 0
    assert summary["exact_gromacs_audit"]["zero_pairs_in_patched_pairs"] == 0
    assert abs(summary["charge"]["charge_mismatch_e"]) <= 1.0e-4
    assert summary["energy"]["openmm_short_dynamics_steps"] == 10
    assert summary["artifacts"]["summary"] == str(summary_path)
    assert Path(summary["artifacts"]["exact_gromacs_audit"]).exists()


def test_final_e2e_summary_validation_rejects_lossy_gromacs_audit() -> None:
    """Summary validation should reject nonzero exact GROMACS mismatch counters."""

    summary = _minimal_valid_summary()
    summary["exact_gromacs_audit"]["pair_mismatch_count"] = 1

    with pytest.raises(AssertionError, match="pair_mismatch_count"):
        _assert_final_e2e_summary(summary)


def test_final_e2e_gromacs_runner_sets_cwd(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """GROMACS subprocess helper should run inside the per-test work directory."""

    calls = []

    def fake_run(**kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append(kwargs)
        return subprocess.CompletedProcess(kwargs["args"], 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    _run(["gmx", "--version"], cwd=tmp_path)

    assert calls[0]["cwd"] == tmp_path


def _require_final_e2e_prerequisites() -> None:
    """Fail clearly when final E2E runtime prerequisites are unavailable."""

    try:
        import openmm  # noqa: F401
    except Exception as exc:  # pragma: no cover - depends on external env
        pytest.fail(f"Final E2E requires OpenMM in the active pixi environment: {exc}")
    if shutil.which("gmx") is None:
        pytest.fail("Final E2E requires a GROMACS 'gmx' executable on PATH")


def _write_final_runtime_config(tmp_path: Path) -> Path:
    """Write a runtime config with the full final component mix."""

    payload = yaml.safe_load(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    payload["name"] = "final_e2e_1ubq"
    payload["enzyme"]["pdb_path"] = str(FIXTURE_DIR / "1ubq.pdb")
    payload["polymers"]["cache_directory"] = str(tmp_path / "polymer_cache")
    payload["polymers"]["count"] = 1
    payload["output"]["projects_directory"] = str(tmp_path / "projects")
    payload["output"]["scratch_directory"] = str(tmp_path / "scratch")
    payload["solvent"]["co_solvents"] = [{"name": "ethanol", "mole_fraction": 0.02}]
    payload["solvent"]["ions"]["nacl_concentration"] = 0.01
    glycan_one = _runtime_glycan_fragment(tmp_path, "glycan_asn60.pdb")
    glycan_two = _runtime_glycan_fragment(tmp_path, "glycan_asn25.pdb")
    attachments = payload["conjugation"]["attachments"]
    attachments[1]["moiety"] = {
        "name": "glygen_g42666_asn60",
        "force_field": "glycam06",
        "input_path": str(glycan_one),
    }
    polymer = yaml.safe_load(yaml.safe_dump(attachments[0]))
    polymer["name"] = "abc_lys6_polymer"
    polymer["site"]["residue_number"] = 6
    glycan = yaml.safe_load(yaml.safe_dump(attachments[1]))
    glycan["name"] = "asn25_three_branch_nglycan"
    glycan["site"]["residue_number"] = 25
    glycan["moiety"] = {
        "name": "glygen_g42666_asn25",
        "force_field": "glycam06",
        "input_path": str(glycan_two),
    }
    attachments.extend([polymer, glycan])
    config_path = tmp_path / "final_e2e_config.yaml"
    config_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return config_path


def _runtime_glycan_fragment(tmp_path: Path, name: str) -> Path:
    """Write a TER-terminated runtime copy of the compact GLYCAM PDB fragment."""

    source = FIXTURE_DIR.parent / "glygen" / "synthetic_g42666_style_conect.pdb"
    text = source.read_text(encoding="utf-8")
    lines = text.splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    conect_lines = [line for line in lines if line.startswith("CONECT")]
    destination = tmp_path / "runtime_inputs" / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        "\n".join([*atom_lines, "TER", *conect_lines, "END", ""]), encoding="utf-8"
    )
    return destination


def _run_short_openmm_exact_bundle(bundle: Any, output_dir: Path) -> dict[str, Any]:
    """Run uncapped minimization and short dynamics with the exact OpenMM bundle."""

    import openmm
    from openmm import unit
    from openmm.app import Simulation

    output_dir.mkdir(parents=True, exist_ok=True)
    system = bundle.to_openmm()
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtosecond,
    )
    simulation = Simulation(bundle.to_openmm_topology(), system, integrator)
    simulation.context.setPositions(bundle.to_openmm_positions())
    simulation.minimizeEnergy(maxIterations=0)
    simulation.step(10)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return {"openmm_short_dynamics_steps": 10, "openmm_potential_energy_kj_mol": float(energy)}


def _run_gromacs_acceptance(exported: dict[str, Any], work_dir: Path) -> dict[str, Any]:
    """Run GROMACS grompp/minimize/rerun acceptance commands."""

    work_dir.mkdir(parents=True, exist_ok=True)
    gmx = shutil.which("gmx") or "gmx"
    top = Path(exported["top"])
    gro = Path(exported["gro"])
    em_mdp = Path(exported["em_mdp"])
    prod_mdp = Path(exported["prod_mdp"])
    tpr = work_dir / "final_e2e_em.tpr"
    prod_tpr = work_dir / "final_e2e_prod.tpr"
    _run(
        [
            gmx,
            "grompp",
            "-f",
            str(em_mdp),
            "-c",
            str(gro),
            "-p",
            str(top),
            "-o",
            str(tpr),
            "-po",
            str(work_dir / "final_e2e_em_out.mdp"),
            "-maxwarn",
            "0",
        ],
        cwd=work_dir,
    )
    _run([gmx, "mdrun", "-deffnm", str(work_dir / "final_e2e_em")], cwd=work_dir)
    _run(
        [
            gmx,
            "grompp",
            "-f",
            str(prod_mdp),
            "-c",
            str(work_dir / "final_e2e_em.gro"),
            "-p",
            str(top),
            "-o",
            str(prod_tpr),
            "-po",
            str(work_dir / "final_e2e_prod_out.mdp"),
            "-maxwarn",
            "0",
        ],
        cwd=work_dir,
    )
    rerun_log = work_dir / "final_e2e_rerun.log"
    _run(
        [
            gmx,
            "mdrun",
            "-s",
            str(prod_tpr),
            "-rerun",
            str(work_dir / "final_e2e_em.trr"),
            "-g",
            str(rerun_log),
            "-e",
            str(work_dir / "final_e2e_rerun.edr"),
        ],
        cwd=work_dir,
    )
    return {
        "gromacs_tpr": str(tpr),
        "gromacs_prod_tpr": str(prod_tpr),
        "gromacs_minimized_gro": str(work_dir / "final_e2e_em.gro"),
        "gromacs_rerun_log": str(rerun_log),
    }


def _run(command: list[str], *, cwd: Path) -> None:
    """Run a subprocess and raise an assertion with captured diagnostics on failure."""

    try:
        completed = subprocess.run(
            args=command,
            text=True,
            capture_output=True,
            check=False,
            timeout=GROMACS_E2E_TIMEOUT_SECONDS,
            cwd=cwd,
        )
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(
            "Command timed out after "
            f"{GROMACS_E2E_TIMEOUT_SECONDS} s in {cwd}: {' '.join(command)}"
            f"\nstdout:\n{exc.stdout or ''}\nstderr:\n{exc.stderr or ''}"
        ) from exc
    if completed.returncode != 0:
        raise AssertionError(
            "Command failed: "
            + " ".join(command)
            + f"\ncwd: {cwd}"
            + f"\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )


def _final_e2e_summary(
    result: Any,
    exported: dict[str, Any],
    openmm_summary: dict[str, Any],
    gromacs_summary: dict[str, Any],
) -> dict[str, Any]:
    """Return final E2E evidence summary."""

    overlay = _load_json(result.artifact_paths["overlay_diagnostics"])
    ownership = _load_json(result.artifact_paths["ownership_manifest"])
    charge = _load_json(result.artifact_paths["mixed_overlay_charge_audit"])
    native_audit = _load_json(result.artifact_paths["native_openmm_glycam_audit"])
    exact_gromacs_audit = _load_json(exported["exact_gromacs_audit"])
    validation = _load_json(result.artifact_paths["conjugate_validation_report"])
    atoms = _pdb_atoms(Path(result.solvated_pdb_path))
    artifacts = {name: str(path) for name, path in result.artifact_paths.items()}
    artifacts.update({name: str(path) for name, path in exported.items() if isinstance(path, Path)})
    artifacts.update(gromacs_summary)
    linkage_evidence = _derive_linkage_evidence(validation)
    component_counts = _derive_component_counts(
        atoms=atoms,
        linkage_evidence=linkage_evidence,
        native_audit=native_audit,
    )
    return {
        "components": component_counts,
        "linkages": linkage_evidence,
        "plan_counts": {
            "attachment_specs": len(result.attachment_specs),
            "glycan_specs": sum(
                1
                for spec in result.attachment_specs
                if getattr(spec, "reaction_name", "") == "n_glycosylation"
                or "pdb" in spec.fragment.sidecars
            ),
        },
        "ownership": {
            "glycam_atoms": len(ownership["domains"]["glycam"]),
            "generic_atoms": len(ownership["domains"]["generic"]),
        },
        "parity": overlay["parity"],
        "exact_gromacs_audit": exact_gromacs_audit,
        "charge": charge,
        "energy": openmm_summary,
        "artifacts": artifacts,
    }


def _load_json(path: Path | str) -> dict[str, Any]:
    """Load a JSON object from a path."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


def _pdb_atoms(path: Path) -> list[dict[str, Any]]:
    """Return compact atom records from a PDB file."""

    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(
                {
                    "serial": int(line[6:11]),
                    "name": line[12:16].strip(),
                    "chain": line[21].strip(),
                    "residue": line[17:20].strip(),
                    "residue_number": int(line[22:26]),
                    "insertion_code": line[26].strip(),
                }
            )
    return atoms


def _derive_linkage_evidence(validation: dict[str, Any]) -> dict[str, Any]:
    """Derive final linkage evidence from the validation bond graph artifact."""

    product_path = Path(validation["product_pdb_path"])
    atoms_by_serial = {atom["serial"]: atom for atom in _pdb_atoms(product_path)}
    expected_pairs = _normalized_pair_set(validation["bond_graph"]["expected_bonds"])
    observed_pairs = _normalized_pair_set(validation["bond_graph"]["observed_bonds"])
    missing_pairs = _normalized_pair_set(validation["bond_graph"].get("missing_bonds", ()))
    present_pairs = tuple(sorted(expected_pairs & observed_pairs))
    classified = [_classify_linkage_pair(pair, atoms_by_serial) for pair in present_pairs]
    glycan_pairs = [entry for entry in classified if entry["classification"] == "glycan"]
    polymer_pairs = [entry for entry in classified if entry["classification"] == "polymer"]
    return {
        "expected_count": len(expected_pairs),
        "observed_expected_count": len(present_pairs),
        "missing_count": len(missing_pairs),
        "glycan_linkages": len(glycan_pairs),
        "polymer_linkages": len(polymer_pairs),
        "pairs": classified,
    }


def _derive_component_counts(
    *,
    atoms: list[dict[str, Any]],
    linkage_evidence: dict[str, Any],
    native_audit: dict[str, Any],
) -> dict[str, int]:
    """Derive final component counts from produced PDB and audit artifacts."""

    sage_components = tuple(native_audit.get("sage_template_generator", {}).get("components", ()))
    free_polymers = _count_sage_components_with_residues(sage_components, POLYMER_RESIDUE_NAMES)
    if free_polymers == 0:
        free_polymers = _count_free_polymer_chains(atoms)
    return {
        "attached_glycans": int(linkage_evidence["glycan_linkages"]),
        "attached_polymers": int(linkage_evidence["polymer_linkages"]),
        "free_polymers": free_polymers,
        "cosolvents": _count_pdb_residues(atoms, COSOLVENT_RESIDUE_NAMES),
        "water_atoms": sum(1 for atom in atoms if atom["residue"] in WATER_RESIDUE_NAMES),
        "ion_atoms": sum(1 for atom in atoms if atom["residue"] in ION_RESIDUE_NAMES),
    }


def _count_free_polymer_chains(atoms: list[dict[str, Any]]) -> int:
    """Count final PDB non-product chains that contain polymer residues."""

    polymer_chains = {
        atom["chain"]
        for atom in atoms
        if atom["residue"] in POLYMER_RESIDUE_NAMES and atom["chain"] not in {"A", "C"}
    }
    return len(polymer_chains)


def _count_sage_components_with_residues(
    components: tuple[dict[str, Any], ...], residue_names: frozenset[str]
) -> int:
    """Count audited Sage components containing at least one target residue name."""

    count = 0
    for component in components:
        segment_names = {
            str(segment.get("residue_name", "")).strip()
            for segment in component.get("source_residue_segments", ())
        }
        if segment_names & residue_names:
            count += 1
    return count


def _count_pdb_residues(atoms: list[dict[str, Any]], residue_names: frozenset[str]) -> int:
    """Count unique final PDB residues matching target residue names."""

    return len(
        {
            (atom["chain"], atom["residue"], atom["residue_number"], atom["insertion_code"])
            for atom in atoms
            if atom["residue"] in residue_names
        }
    )


def _normalized_pair_set(pairs: Any) -> set[tuple[int, int]]:
    """Return sorted integer pair tuples from JSON pair records."""

    return {tuple(sorted((int(pair[0]), int(pair[1])))) for pair in pairs}


def _classify_linkage_pair(
    pair: tuple[int, int], atoms_by_serial: dict[int, dict[str, Any]]
) -> dict[str, Any]:
    """Classify one linkage pair from product PDB atom identities."""

    left = atoms_by_serial[pair[0]]
    right = atoms_by_serial[pair[1]]
    residues = {left["residue"], right["residue"]}
    if residues & GLYCAM_RESIDUE_NAMES:
        classification = "glycan"
    elif residues & POLYMER_RESIDUE_NAMES:
        classification = "polymer"
    else:
        classification = "other"
    return {"pair": pair, "classification": classification, "atoms": (left, right)}


def _assert_final_e2e_summary(summary: dict[str, Any]) -> None:
    """Validate final E2E summary evidence before and after serialization."""

    assert summary["linkages"]["expected_count"] == 4, "expected_count"
    assert summary["linkages"]["observed_expected_count"] == 4, "observed_expected_count"
    assert summary["linkages"]["missing_count"] == 0, "missing_count"
    assert summary["linkages"]["glycan_linkages"] == 2, "glycan_linkages"
    assert summary["linkages"]["polymer_linkages"] == 2, "polymer_linkages"
    assert summary["exact_gromacs_audit"]["pair_mismatch_count"] == 0, "pair_mismatch_count"
    assert (
        summary["exact_gromacs_audit"]["exclusion_mismatch_count"] == 0
    ), "exclusion_mismatch_count"
    assert (
        summary["exact_gromacs_audit"]["zero_pairs_in_patched_pairs"] == 0
    ), "zero_pairs_in_patched_pairs"


def _minimal_valid_summary() -> dict[str, Any]:
    """Return a compact summary object for validation unit tests."""

    return {
        "linkages": {
            "expected_count": 4,
            "observed_expected_count": 4,
            "missing_count": 0,
            "glycan_linkages": 2,
            "polymer_linkages": 2,
        },
        "exact_gromacs_audit": {
            "pair_mismatch_count": 0,
            "exclusion_mismatch_count": 0,
            "zero_pairs_in_patched_pairs": 0,
        },
    }
