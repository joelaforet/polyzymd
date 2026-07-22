"""Opt-in fast mixed polymer/glycan conjugation E2E helper."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any, Literal

import pytest
import yaml

from polyzymd.builders.conjugation.api import build_conjugate_from_config
from polyzymd.builders.conjugation.pablo.parameterization import InterchangeParameterizationSettings
from polyzymd.builders.conjugation.relaxation import ConjugateRelaxationSettings
from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings
from polyzymd.config.loader import load_config
from polyzymd.exporters.gromacs import GromacsExporter

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "data/conjugation/fast_mixed_e2e_1ubq"
FIXTURE_CONFIG = FIXTURE_DIR / "config.yaml"
SUMMARY_NAME = "fast_mixed_e2e_summary.json"
GlycosylationMode = Literal["n", "o_ser", "o_thr"]


def run_fast_mixed_build_export(
    output_dir: Path,
    *,
    run_cuda: bool = False,
    cuda_steps: int = 100_000,
    glycosylation: GlycosylationMode = "n",
) -> dict[str, Any]:
    """Run the fast mixed conjugation checkpoint and return a compact summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    runtime_config = _runtime_config_path(output_dir, glycosylation=glycosylation)
    config = load_config(runtime_config)
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
    start = time.perf_counter()
    result = build_conjugate_from_config(
        config,
        output_dir=output_dir / "build",
        settings=settings,
        free_polymer_seed=202610,
    )
    parameterized_system = result.exact_export_bundle or result.final_interchange
    if parameterized_system is None:
        raise AssertionError("Fast mixed E2E build did not create a parameterized system")
    exported = GromacsExporter(parameterized_system, config).export(
        output_dir / "gromacs",
        prefix="fast_mixed_e2e_1ubq",
    )
    summary = _build_summary(
        result,
        exported=exported,
        build_seconds=time.perf_counter() - start,
        glycosylation=glycosylation,
    )
    if run_cuda:
        summary["cuda"] = _run_cuda_smoke(parameterized_system, output_dir / "cuda", cuda_steps)
    (output_dir / SUMMARY_NAME).write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def validate_fast_mixed_summary(summary: dict[str, Any], *, require_cuda: bool = False) -> None:
    """Assert that the compact summary satisfies the checkpoint contract."""
    assert summary["fixture_glycan"]["heavy_atom_count"] == 117
    assert summary["fixture_glycan"]["hydrogenated_atom_count"] == 225
    assert summary["fixture_glycan"]["formula"] == "C64H108N2O51"
    assert summary["fixture_glycan"]["anomeric_carbon_index"] == 9
    assert summary["fixture_glycan"]["leaving_atom_indices"] == [10, 126]
    assert summary["fixture_glycan"]["leaving_atom_count"] == 2
    assert summary["authoritative_system_created"] is True
    assert summary["resolved_attachment_count"] == 2
    assert summary["solvated_atom_count"] > summary["crosslinked_atom_count"] > 1_000
    assert summary["finite_coordinates"] is True
    assert summary["residue_counts"]["LYX"] > 0
    assert summary["residue_counts"]["ASX"] > 0
    assert summary["residue_counts"]["NAG"] == 223
    assert summary["product_glycan"]["residue_name"] == "NAG"
    assert summary["product_glycan"]["source_atom_count"] == 225
    assert summary["product_glycan"]["leaving_atom_count"] == 2
    assert summary["product_glycan"]["retained_atom_count"] == 223
    assert summary["product_glycan"]["retained_atom_count"] == (
        summary["product_glycan"]["source_atom_count"]
        - summary["product_glycan"]["leaving_atom_count"]
    )
    n_glycan = summary["n_glycosylation_linkage"]
    assert n_glycan["mechanism_name"] == "n_glycosylation"
    assert n_glycan["site_residue_number"] == 60
    assert n_glycan["protein_product_residue_name"] == "ASX"
    assert n_glycan["modifier_product_residue_name"] == "NAG"
    assert n_glycan["protein_link_atom_name"] == "ND2"
    assert n_glycan["product_site_atom_present"] is True
    assert n_glycan["protein_leaving_atom_names"] == ["HD21"]
    assert n_glycan["modifier_leaving_atom_count"] == 2
    assert n_glycan["crosslink_residues"] == ["ASX", "NAG"]
    assert n_glycan["crosslink_linking_atoms"][0] == "ND2"
    assert all(summary["export_exists"].values())
    validation = summary["validation"]
    assert validation["bond_graph_status"] == "pass"
    assert validation["atom_presence_status"] == "pass"
    assert validation["charge_audit_status"] == "pass"
    assert validation["relaxation_evidence_status"] == "pass"
    assert len(validation["linkage_distances_angstrom"]) == 2
    assert validation["close_contact_count"] == 0
    charge_bridge = summary["charge_bridge"]
    assert charge_bridge["success"] is True
    assert charge_bridge["ff14sb_atom_count"] > 1_000
    assert charge_bridge["local_nagl_patch_atom_count"] > 0
    assert charge_bridge["source_role_atom_count"] == summary["crosslinked_atom_count"]
    assert charge_bridge["polymer_atom_count"] > 0
    assert charge_bridge["polymer_source_coverage_count"] == charge_bridge["polymer_atom_count"]
    assert charge_bridge["patch_owned_polymer_atom_count"] == charge_bridge["polymer_atom_count"]
    assert (
        charge_bridge["polymer_template_coverage_count"]
        == charge_bridge["polymer_template_mapped_product_atom_count"]
    )
    assert (
        charge_bridge["polymer_template_atom_count"]
        == charge_bridge["unmodified_polymer_template_atom_count"]
    )
    assert charge_bridge["polymer_template_overwrite_count"] > 0
    assert charge_bridge["invalid_polymer_template_source_count"] == 0
    assert charge_bridge["total_charge_e"] == pytest.approx(charge_bridge["formal_charge_e"])
    assert abs(charge_bridge["normalization_correction_e"] or 0.0) <= 0.02
    assert abs(charge_bridge["max_per_atom_correction_e"] or 0.0) <= 0.005
    relaxation = summary["relaxation"]
    assert math.isfinite(relaxation["stage_a_energy_before_min_kj_mol"])
    assert math.isfinite(relaxation["stage_a_energy_after_min_kj_mol"])
    if require_cuda:
        cuda = summary["cuda"]
        assert cuda["platform"] == "CUDA"
        assert cuda["steps"] == 100_000
        assert cuda["checkpoint_exists"] is True
        assert math.isfinite(cuda["potential_energy_kj_mol"])


def validate_fast_mixed_o_summary(
    summary: dict[str, Any],
    *,
    product_residue: Literal["OLS", "OLT"] = "OLS",
    residue_number: int = 60,
    site_atom: Literal["OG", "OG1"] = "OG",
    require_cuda: bool = False,
) -> None:
    """Assert a Ser or Thr O-glycosylation mixed checkpoint contract."""
    assert (product_residue, site_atom) in {("OLS", "OG"), ("OLT", "OG1")}
    assert summary["fixture_glycan"]["source"] == "glycam_G57321FI_CONECT.pdb"
    assert summary["fixture_glycan"]["hydrogenated_atom_count"] == 30
    assert summary["fixture_glycan"]["anomeric_carbon_serial"] == 4
    assert summary["fixture_glycan"]["leaving_atom_serials"] == [1, 2]
    assert summary["authoritative_system_created"] is True
    assert summary["resolved_attachment_count"] == 2
    assert summary["solvated_atom_count"] > summary["crosslinked_atom_count"] > 1_000
    assert summary["finite_coordinates"] is True
    assert summary["residue_counts"]["LYX"] > 0
    assert summary["residue_counts"][product_residue] > 0
    assert summary["residue_counts"]["0VA"] == 28
    linkage = summary["glycosylation_linkage"]
    assert linkage["mechanism_name"] == "o_glycosylation"
    assert linkage["site_residue_number"] == residue_number
    assert linkage["protein_product_residue_name"] == product_residue
    assert linkage["modifier_product_residue_name"] == "0VA"
    assert linkage["protein_link_atom_name"] == site_atom
    assert linkage["product_site_atom_present"] is True
    expected_leaving_hydrogen = "HG" if site_atom == "OG" else "HG1"
    assert linkage["protein_leaving_atom_names"] == [expected_leaving_hydrogen]
    assert linkage["modifier_leaving_atom_count"] == 2
    assert linkage["crosslink_residues"] == [product_residue, "0VA"]
    assert linkage["crosslink_linking_atoms"][0] == site_atom
    assert all(summary["export_exists"].values())
    validation = summary["validation"]
    assert validation["bond_graph_status"] == "pass"
    assert validation["atom_presence_status"] == "pass"
    assert validation["charge_audit_status"] == "pass"
    assert validation["relaxation_evidence_status"] == "pass"
    assert len(validation["linkage_distances_angstrom"]) == 2
    assert validation["close_contact_count"] == 0
    assert summary["charge_bridge"]["success"] is True
    if require_cuda:
        cuda = summary["cuda"]
        assert cuda["platform"] == "CUDA"
        assert cuda["steps"] > 0
        assert cuda["checkpoint_exists"] is True
        assert math.isfinite(cuda["potential_energy_kj_mol"])


def _runtime_config_path(
    output_dir: Path,
    *,
    glycosylation: GlycosylationMode = "n",
) -> Path:
    """Write a runtime-local config so fixture directories stay read-only."""
    payload = yaml.safe_load(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    runtime_dir = output_dir / "runtime_inputs"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    protein_path = FIXTURE_DIR / "1ubq.pdb"
    if glycosylation in {"o_ser", "o_thr"}:
        protein_path = runtime_dir / f"1ubq_{glycosylation}60.pdb"
        protein_builder = (
            _ser60_runtime_protein if glycosylation == "o_ser" else _thr60_runtime_protein
        )
        protein_path.write_text(
            protein_builder(FIXTURE_DIR / "1ubq.pdb"),
            encoding="utf-8",
        )
    if glycosylation != "n":
        glycan_attachment = payload["conjugation"]["attachments"][1]
        if glycosylation == "o_ser":
            glycan_attachment["name"] = "ser60_g57321fi_oglycan"
            glycan_attachment["site"].update({"residue_name": "SER", "atom_name": "OG"})
        else:
            glycan_attachment["name"] = "thr60_g57321fi_oglycan"
            glycan_attachment["site"].update(
                {"residue_name": "THR", "residue_number": 60, "atom_name": "OG1"}
            )
        glycan_attachment["moiety"] = {
            "name": "glycam_G57321FI",
            "force_field": "glycam06",
            "input_path": str(FIXTURE_DIR / "glycam_G57321FI_CONECT.pdb"),
        }
        glycan_attachment["mechanism"]["name"] = "o_glycosylation"
    payload["enzyme"]["pdb_path"] = str(protein_path)
    payload["polymers"]["cache_directory"] = str(runtime_dir / "polymer_cache")
    payload["output"]["projects_directory"] = str(output_dir / "projects")
    payload["output"]["scratch_directory"] = str(output_dir / "scratch")
    runtime_config = runtime_dir / "config.yaml"
    runtime_config.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return runtime_config


def _build_summary(
    result: Any,
    *,
    exported: dict[str, Any],
    build_seconds: float,
    glycosylation: GlycosylationMode = "n",
) -> dict[str, Any]:
    """Return a compact JSON-compatible build/export validation summary."""
    solvated_pdb = Path(result.solvated_pdb_path)
    crosslinked_pdb = Path(result.construction.crosslinked_pdb_path)
    validation = _load_json(getattr(result.construction, "validation_report_path", None))
    charge_bridge = _load_json(
        Path(result.construction.output_dir) / "product_state_charge_bridge.json"
    )
    relaxation = _load_json(Path(result.construction.output_dir) / "conjugate_relaxation.json")
    atoms = _pdb_atom_records(solvated_pdb)
    crosslinked_atoms = _pdb_atom_records(crosslinked_pdb)
    export_paths = _flatten_paths(exported)
    resolved = [product.model_dump(mode="json") for product in result.attachment_specs]
    mechanism_name = "n_glycosylation" if glycosylation == "n" else "o_glycosylation"
    product_site = {"n": "ASX", "o_ser": "OLS", "o_thr": "OLT"}[glycosylation]
    site_atom = {"n": "ND2", "o_ser": "OG", "o_thr": "OG1"}[glycosylation]
    site_residue_number = 60
    summary = {
        "fixture_config": str(FIXTURE_CONFIG),
        "fixture_glycan": _glycan_fixture_evidence(glycosylation=glycosylation),
        "build_seconds": build_seconds,
        "solvated_pdb": str(solvated_pdb),
        "crosslinked_pdb": str(crosslinked_pdb),
        "solvated_atom_count": len(atoms),
        "crosslinked_atom_count": len(crosslinked_atoms),
        "finite_coordinates": all(
            math.isfinite(coord) for atom in atoms for coord in (atom["x"], atom["y"], atom["z"])
        ),
        "authoritative_system_created": bool(
            result.exact_export_bundle is not None or result.final_interchange_created
        ),
        "resolved_attachment_count": len(resolved),
        "resolved_plans": resolved,
        "residue_counts": _residue_counts(crosslinked_atoms),
        "product_glycan": _product_glycan_evidence(crosslinked_atoms, glycosylation=glycosylation),
        "glycosylation_linkage": _glycosylation_linkage_evidence(
            resolved,
            crosslinked_atoms,
            mechanism_name=mechanism_name,
            product_site_residue_name=product_site,
            site_residue_number=site_residue_number,
            site_atom_name=site_atom,
        ),
        "validation": _validation_summary(validation),
        "charge_bridge": _charge_bridge_summary(charge_bridge, crosslinked_atoms),
        "relaxation": _relaxation_summary(relaxation),
        "exports": {name: str(path) for name, path in export_paths.items()},
        "export_exists": {name: path.exists() for name, path in export_paths.items()},
    }
    if glycosylation == "n":
        summary["n_glycosylation_linkage"] = summary["glycosylation_linkage"]
    return summary


def _glycan_fixture_evidence(*, glycosylation: GlycosylationMode = "n") -> dict[str, Any]:
    """Return source evidence for the selected glycan fixture."""
    if glycosylation != "n":
        atoms = _pdb_atom_records(FIXTURE_DIR / "glycam_G57321FI_CONECT.pdb")
        return {
            "source": "glycam_G57321FI_CONECT.pdb",
            "hydrogenated_atom_count": len(atoms),
            "anomeric_carbon_serial": 4,
            "leaving_atom_serials": [1, 2],
            "leaving_atom_count": 2,
        }
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    from polyzymd.builders.conjugation.reactions.n_glycosylation import (
        detect_glycan_anomeric_group,
    )

    payload = yaml.safe_load(FIXTURE_CONFIG.read_text(encoding="utf-8"))
    smiles = payload["conjugation"]["attachments"][1]["moiety"]["smiles"]
    heavy_mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(heavy_mol)
    group = detect_glycan_anomeric_group(mol)
    return {
        "heavy_atom_count": heavy_mol.GetNumAtoms(),
        "hydrogenated_atom_count": mol.GetNumAtoms(),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "anomeric_carbon_index": group.reactive_carbon_index,
        "leaving_atom_indices": list(group.leaving_atom_indices),
        "leaving_atom_count": len(group.leaving_atom_indices),
    }


def _product_glycan_evidence(
    atoms: list[dict[str, Any]], *, glycosylation: GlycosylationMode = "n"
) -> dict[str, Any]:
    """Return retained product glycan atom-count evidence from the product PDB."""
    fixture = _glycan_fixture_evidence(glycosylation=glycosylation)
    residue_name = "NAG" if glycosylation == "n" else "0VA"
    retained_atom_count = sum(1 for atom in atoms if atom["residue_name"] == residue_name)
    leaving_atom_count = int(fixture["leaving_atom_count"])
    return {
        "residue_name": residue_name,
        "source_atom_count": fixture["hydrogenated_atom_count"],
        "leaving_atom_count": leaving_atom_count,
        "retained_atom_count": retained_atom_count,
    }


def _glycosylation_linkage_evidence(
    resolved_plans: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    *,
    mechanism_name: str,
    product_site_residue_name: str,
    site_residue_number: int,
    site_atom_name: str,
) -> dict[str, Any]:
    """Return resolved N- or O-glycosylation linkage evidence."""
    plan = _find_resolved_plan(resolved_plans, mechanism_name)
    protein_link = plan["protein_link_atom"]
    requirement = plan["pablo_crosslink_requirement"]
    return {
        "mechanism_name": plan["contract"]["mechanism_name"],
        "site_residue_number": protein_link["residue_number"],
        "protein_product_residue_name": plan["protein_product_residue_name"],
        "modifier_product_residue_name": plan["modifier_product_residue_name"],
        "protein_link_atom_name": protein_link["atom_name"],
        "product_site_atom_present": _has_product_atom(
            atoms, product_site_residue_name, site_residue_number, site_atom_name
        ),
        "protein_leaving_atom_names": [
            atom["atom_name"] for atom in plan.get("protein_leaving_atoms", [])
        ],
        "modifier_leaving_atom_count": len(plan.get("modifier_leaving_atoms", [])),
        "crosslink_residues": list(requirement["residues"]),
        "crosslink_linking_atoms": list(requirement["linking_atoms"]),
    }


def _ser60_runtime_protein(source_path: Path) -> str:
    """Return the 1UBQ fixture with ASN60 replaced by a standard Ser sidechain."""
    output = []
    for line in source_path.read_text(encoding="utf-8").splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            output.append(line)
            continue
        if line[21:22] != "A" or int(line[22:26]) != 60:
            output.append(line)
            continue
        atom_name = line[12:16].strip()
        if atom_name in {"OD1", "ND2"}:
            continue
        updated = f"{line[:17]}SER{line[20:]}"
        if atom_name == "CG":
            updated = f"{updated[:12]} OG {updated[16:76]} O  {updated[80:]}"
        output.append(updated)
    return "".join(output)


def _thr60_runtime_protein(source_path: Path) -> str:
    """Return 1UBQ with ASN60 replaced by a complete heavy-atom Thr sidechain."""
    lines = source_path.read_text(encoding="utf-8").splitlines(keepends=True)
    residue_lines = [
        line
        for line in lines
        if line.startswith(("ATOM  ", "HETATM")) and line[21:22] == "A" and int(line[22:26]) == 60
    ]
    coordinates = {
        line[12:16].strip(): tuple(
            float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54))
        )
        for line in residue_lines
    }
    cb = coordinates["CB"]
    bond_vectors = []
    for neighbor_name in ("CA", "CG"):
        vector = tuple(coordinates[neighbor_name][index] - cb[index] for index in range(3))
        vector_norm = math.sqrt(sum(value * value for value in vector))
        bond_vectors.append(tuple(value / vector_norm for value in vector))
    direction = tuple(-(bond_vectors[0][index] + bond_vectors[1][index]) for index in range(3))
    norm = math.sqrt(sum(value * value for value in direction))
    cg2 = tuple(cb[index] + 1.53 * direction[index] / norm for index in range(3))
    output = []
    for line in lines:
        if (
            not line.startswith(("ATOM  ", "HETATM"))
            or line[21:22] != "A"
            or int(line[22:26]) != 60
        ):
            output.append(line)
            continue
        atom_name = line[12:16].strip()
        if atom_name == "ND2":
            continue
        updated = f"{line[:17]}THR{line[20:]}"
        if atom_name == "CG":
            updated = f"{updated[:12]} OG1{updated[16:76]} O  {updated[80:]}"
        elif atom_name == "OD1":
            updated = (
                f"{updated[:12]} CG2{updated[16:30]}"
                f"{cg2[0]:8.3f}{cg2[1]:8.3f}{cg2[2]:8.3f}{updated[54:76]} C  {updated[80:]}"
            )
        output.append(updated)
    return "".join(output)


def _find_resolved_plan(
    resolved_plans: list[dict[str, Any]], mechanism_name: str
) -> dict[str, Any]:
    """Find the resolved attachment plan for a named reaction mechanism."""
    for plan in resolved_plans:
        if plan.get("contract", {}).get("mechanism_name") == mechanism_name:
            return plan
    raise AssertionError(f"No resolved attachment plan for mechanism {mechanism_name!r}")


def _has_product_atom(
    atoms: list[dict[str, Any]], residue_name: str, residue_number: int, atom_name: str
) -> bool:
    """Return whether the product PDB contains a specific residue atom."""
    return any(
        atom["residue_name"] == residue_name
        and atom["residue_number"] == residue_number
        and atom["atom_name"] == atom_name
        for atom in atoms
    )


def _run_cuda_smoke(parameterized_system: Any, output_dir: Path, steps: int) -> dict[str, Any]:
    """Run the optional CUDA smoke test without platform fallback."""
    import openmm
    from openmm import unit
    from openmm.app import CheckpointReporter, Simulation, StateDataReporter

    output_dir.mkdir(parents=True, exist_ok=True)
    platform = openmm.Platform.getPlatformByName("CUDA")
    if getattr(parameterized_system, "is_exact_export_bundle", False):
        system = parameterized_system.to_openmm()
        topology = parameterized_system.to_openmm_topology()
        positions = parameterized_system.to_openmm_positions()
    else:
        system = parameterized_system.to_openmm(combine_nonbonded_forces=True)
        topology = parameterized_system.to_openmm_topology()
        positions = parameterized_system.positions.to_openmm()
    simulation = Simulation(
        topology,
        system,
        openmm.LangevinMiddleIntegrator(
            300.0 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond
        ),
        platform,
    )
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy(maxIterations=0)
    checkpoint_path = output_dir / "cuda.chk"
    simulation.reporters.append(CheckpointReporter(str(checkpoint_path), max(1, steps // 2)))
    simulation.reporters.append(
        StateDataReporter(
            str(output_dir / "cuda_state.csv"), max(1, steps // 10), step=True, potentialEnergy=True
        )
    )
    simulation.step(steps)
    energy = (
        simulation.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )
    return {
        "platform": platform.getName(),
        "steps": steps,
        "potential_energy_kj_mol": float(energy),
        "checkpoint": str(checkpoint_path),
        "checkpoint_exists": checkpoint_path.exists(),
    }


def _validation_summary(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract lightweight validation report fields for assertions."""
    linkage = payload.get("linkage_geometry", {}) or {}
    return {
        "status": payload.get("status"),
        "bond_graph_status": (payload.get("bond_graph", {}) or {}).get("status"),
        "atom_presence_status": (payload.get("atom_presence", {}) or {}).get("status"),
        "charge_audit_status": (payload.get("charge_audit", {}) or {}).get("status"),
        "relaxation_evidence_status": (payload.get("relaxation_evidence", {}) or {}).get("status"),
        "linkage_distances_angstrom": linkage.get("linkage_distances_angstrom", []),
        "close_contact_count": linkage.get("close_contact_count"),
    }


def _charge_bridge_summary(
    payload: dict[str, Any], crosslinked_atoms: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    """Extract lightweight charge bridge fields for assertions."""
    summary = {
        key: payload.get(key)
        for key in (
            "success",
            "ff14sb_atom_count",
            "polymer_template_atom_count",
            "local_nagl_patch_atom_count",
            "total_charge_e",
            "formal_charge_e",
            "normalization_correction_e",
            "max_per_atom_correction_e",
        )
    }
    diagnostics = payload.get("diagnostic_details", {}) or {}
    patch_ledgers = diagnostics.get("patch_ledgers", []) or []
    patch_overwrites = diagnostics.get("patch_overwrites", []) or []
    polymer_ledgers = diagnostics.get("polymer_ledgers", []) or []
    polymer_template_count = int(payload.get("polymer_template_atom_count") or 0)
    polymer_template_overwrite_count = sum(
        1 for entry in patch_overwrites if entry.get("old_source_role") == "polymer_template"
    )
    mapped_template_count = sum(
        int(entry.get("mapped_product_atom_count") or 0) for entry in polymer_ledgers
    )
    patch_owned_polymer_identities = {
        identity
        for ledger in patch_ledgers
        for identity in ledger.get("affected_atom_identities", []) or []
        if str(identity).startswith("chain C ")
    }
    polymer_atom_count = sum(1 for atom in crosslinked_atoms or [] if atom.get("chain_id") == "C")
    source_role_atom_count = sum(
        int(payload.get(key) or 0)
        for key in (
            "ff14sb_atom_count",
            "polymer_template_atom_count",
            "local_nagl_patch_atom_count",
        )
    )
    summary.update(
        {
            "source_role_atom_count": source_role_atom_count,
            "polymer_atom_count": polymer_atom_count,
            "patch_owned_polymer_atom_count": len(patch_owned_polymer_identities),
            "polymer_source_coverage_count": len(patch_owned_polymer_identities)
            + polymer_template_count,
            "polymer_template_mapped_product_atom_count": mapped_template_count,
            "polymer_template_overwrite_count": polymer_template_overwrite_count,
            "polymer_template_coverage_count": polymer_template_overwrite_count
            + polymer_template_count,
            "unmodified_polymer_template_atom_count": max(
                0,
                mapped_template_count - polymer_template_overwrite_count,
            ),
            "invalid_polymer_template_source_count": _invalid_polymer_template_source_count(
                polymer_ledgers
            ),
        }
    )
    return summary


def _invalid_polymer_template_source_count(polymer_ledgers: list[dict[str, Any]]) -> int:
    """Return polymer template ledgers without charged SDF provenance."""
    return sum(1 for ledger in polymer_ledgers if not ledger.get("charged_sdf"))


def _relaxation_summary(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract Stage A and energy fields from relaxation diagnostics."""
    return {key: payload.get(key) for key in payload if key.startswith("stage_a_")}


def _pdb_atom_records(path: Path) -> list[dict[str, Any]]:
    """Parse ATOM and HETATM records needed by E2E assertions."""
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            records.append(
                {
                    "atom_name": line[12:16].strip(),
                    "residue_name": line[17:20].strip(),
                    "chain_id": line[21].strip(),
                    "residue_number": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            )
    return records


def _residue_counts(atoms: list[dict[str, Any]]) -> dict[str, int]:
    """Count atom records by residue name."""
    counts: dict[str, int] = {}
    for atom in atoms:
        counts[str(atom["residue_name"])] = counts.get(str(atom["residue_name"]), 0) + 1
    return counts


def _flatten_paths(exported: dict[str, Any]) -> dict[str, Path]:
    """Flatten exporter path payloads to runtime files."""
    paths = {
        key: Path(exported[key]) for key in ("gro", "top", "em_mdp", "prod_mdp") if key in exported
    }
    for index, path in enumerate(exported.get("eq_mdps", []) or [], start=1):
        paths[f"eq_mdp_{index}"] = Path(path)
    return paths


def _load_json(path: Path | str | None) -> dict[str, Any]:
    """Load a JSON object if present, otherwise return an empty dict."""
    if path is None or not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))
