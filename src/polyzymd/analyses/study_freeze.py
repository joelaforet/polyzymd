"""Freeze a study folder for publication: ``polyzymd study freeze``.

Freezing never stops for something missing; every gap is a warning. It:

1. checks the publishing metadata (:mod:`~polyzymd.analyses.study_metadata`),
   the git state, and whether each listed analysis's stored results still
   match the study (config hashes, equilibration window, stride, function
   hashes, settings and PolyzyMD version), without loading trajectories;
2. for every replicate whose trajectories are on this machine, hashes its
   trajectory and topology files (SHA-256, cached by path, size and
   modification time), and writes gzipped copies of its engine inputs (the
   OpenMM system XML and topology, or the GROMACS ``.tpr``, ``.top``,
   ``.itp`` and ``.mdp`` files) and its final frame to ``deposit/``;
3. writes ``system_summary.csv`` (box, atoms, waters, ions and composition of
   each replicate), ``manifest.json``, ``md_checklist.yaml`` (the
   Communications Biology reliability and reproducibility checklist, filled
   from the manifest), ``CITATION.cff`` and ``.zenodo.json``;
4. commits those files and the stored results in ``results/`` (outputs of
   PolyzyMD, needed to redraw the figures), never your own uncommitted
   inputs, and tags the commit;
5. lays out ``deposit/``: the tagged study, the engine inputs and final
   frames, with the manifest, README and ``CITATION.cff`` at the top, and
   prepares the upload (:mod:`~polyzymd.analyses.study_upload_guide`):
   ``deposit/upload/`` with exactly the files to add to Zenodo,
   ``deposit/trajectories.csv`` and the step-by-step ``deposit/UPLOAD.md``.
   PolyzyMD uploads and publishes nothing.

``deposit/`` is gitignored. Trajectories are not copied: they are listed in
the manifest by size and SHA-256 and deposited on their own, with their DOIs
in ``metadata.related.trajectories``.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

DEPOSIT = "deposit"
MANIFEST = "manifest.json"
CHECKLIST = "md_checklist.yaml"
SUMMARY = "system_summary.csv"
CITATION = "CITATION.cff"
ZENODO = ".zenodo.json"
#: Files freeze writes in the study folder and commits.
GENERATED = (MANIFEST, CHECKLIST, SUMMARY, CITATION, ZENODO)
MANIFEST_SCHEMA = "polyzymd-study-manifest/1"
#: JSON Schema of the manifest, shipped with PolyzyMD and written into every deposit.
MANIFEST_SCHEMA_FILE = "manifest-1.schema.json"
_IONS = {"NA", "CL", "K", "MG", "ZN", "CA", "SOD", "CLA", "POT", "NA+", "CL-", "K+", "MG2+"}
_HASH_CACHE = ".hashes.json"


@dataclass
class FreezeResult:
    """What :func:`freeze` did."""

    root: Path
    tag: str | None
    commit: str | None
    deposit: Path
    manifest: dict[str, Any]
    warnings: list[str] = field(default_factory=list)
    guide: Path | None = None
    upload: Path | None = None


class _Hashes:
    """SHA-256 of files and their sizes, from the shared hash cache (:mod:`~polyzymd.analyses.shared.file_hashes`)."""

    def __init__(self, path: Path) -> None:
        self.path = path
        try:
            self.cache = json.loads(path.read_text())
        except (OSError, ValueError):
            self.cache = {}

    def __call__(self, file: Path) -> dict[str, Any]:
        from polyzymd.analyses.shared.file_hashes import file_sha256

        # The shared cache, so files hashed by analyses are not hashed again.
        return {"size": file.stat().st_size, "sha256": file_sha256(file)}

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.cache))


def _git(root: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), *arguments], capture_output=True, text=True, timeout=120
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def _versions() -> dict[str, str | None]:
    import platform

    import polyzymd

    versions: dict[str, str | None] = {
        "polyzymd": polyzymd.__version__,
        "python": platform.python_version(),
    }
    for module in (
        "MDAnalysis",
        "numpy",
        "scipy",
        "mdtraj",
        "pymbar",
        "openmm",
        "openff.toolkit",
        "openff.interchange",
    ):
        try:
            imported = __import__(module, fromlist=["__version__"])
            versions[module] = str(
                getattr(imported, "__version__", None) or getattr(imported, "version", None)
            )
        except Exception:  # noqa: BLE001 - an absent or broken package is recorded as absent
            versions[module] = None
    return versions


def _function_hash(record: dict[str, Any], entry: Any) -> str | None:
    """Return the current hash of the function a stored record names, or ``None`` if unknown."""
    import importlib

    from polyzymd.analyses.timeseries import _function_record

    try:
        if entry is not None and entry.function is not None:
            from polyzymd.analyses.user_functions import load_function

            function = load_function(entry.function.file, entry.function.qualname)
        else:
            function = importlib.import_module(record["function"]["module"])
            for part in record["function"]["qualname"].split("."):
                function = getattr(function, part)
        return _function_record(function)["hash"]
    except Exception:  # noqa: BLE001 - an unimportable function is reported, not raised
        return None


def stale_runs(protocol: Any) -> dict[str, list[str]]:
    """Return, for each analysis run, why its stored results may not match the study now.

    A run with no stored results, or whose stored records differ from the
    study in config hash, equilibration window, stride or function hash, or
    whose report differs in a listed setting or in the PolyzyMD version, has
    one reason per difference. No trajectory is read.
    """
    import polyzymd
    from polyzymd.analyses.identity import compute_config_hash
    from polyzymd.config.schema import SimulationConfig

    hashes = {}
    for label, path in protocol.conditions.items():
        try:
            hashes[label] = compute_config_hash(SimulationConfig.from_yaml(path))
        except (OSError, ValueError):
            hashes[label] = None
    reasons: dict[str, list[str]] = {}
    for run, entry in protocol.analyses.items():
        folder = protocol.results_dir(run)
        found: list[str] = []
        records = sorted(folder.glob("polyzymd_results/*/*/replicate_*/record.json"))
        report_path = folder / "report.json"
        if not records or not report_path.is_file():
            reasons[run] = [f"no stored results; run polyzymd analyze {run} --study"]
            continue
        current_hashes: dict[tuple, str | None] = {}
        for path in records:
            record = json.loads(path.read_text())
            label = record.get("condition")
            if label in hashes and hashes[label] and record.get("config_hash") != hashes[label]:
                found.append(f"the config of {label} changed")
            if record.get("equilibration") != protocol.equilibration:
                found.append(
                    f"equilibration {record.get('equilibration')} is not {protocol.equilibration}"
                )
            if int(record.get("stride", 1)) != protocol.stride:
                found.append(f"stride {record.get('stride')} is not {protocol.stride}")
            key = (record["function"]["module"], record["function"]["qualname"])
            if key not in current_hashes:
                current_hashes[key] = _function_hash(record, entry)
            if current_hashes[key] and current_hashes[key] != record["function"]["hash"]:
                found.append(f"the code of {key[1]} changed")
        report = json.loads(report_path.read_text())
        provenance = report.get("provenance", {})
        if provenance.get("polyzymd_version") != polyzymd.__version__:
            found.append(
                f"made with PolyzyMD {provenance.get('polyzymd_version')}, not {polyzymd.__version__}"
            )
        expected = entry.function.settings if entry.function is not None else entry.settings
        ran_with = (provenance.get("study") or {}).get("settings")
        if ran_with is not None:
            for key in sorted(set(expected) | set(ran_with)):
                if json.dumps(expected.get(key), sort_keys=True) != json.dumps(
                    ran_with.get(key), sort_keys=True
                ):
                    found.append(f"setting {key} differs from the one the results were made with")
        if found:
            reasons[run] = sorted(set(found))
    return reasons


def _gzip_copy(source: Path, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as raw, gzip.open(target, "wb") as packed:
        shutil.copyfileobj(raw, packed)
    return target


def _engine_inputs(provenance: Any) -> list[Path]:
    """Return the engine input files of a replicate: OpenMM system XML and topology, or GROMACS files."""
    from polyzymd.analyses.shared.loader import openmm_system_file

    topology = Path(provenance.topology.path)
    files: list[Path] = []
    if (provenance.config_engine or "openmm") == "gromacs":
        folder = topology.parent
        for pattern in ("prod.tpr", "*.top", "*.itp", "*.mdp"):
            files.extend(sorted(folder.glob(pattern)))
    else:
        files.append(topology)
        for trajectory in provenance.trajectories:
            system = openmm_system_file(trajectory.path)
            if system is not None:
                files.append(system)
    build = Path(provenance.working_directory) / "build_manifest.json"
    if build.is_file():
        files.append(build)
    unique, seen = [], set()
    for path in files:
        if path.is_file() and path.resolve() not in seen:
            seen.add(path.resolve())
            unique.append(path)
    return unique


#: What a deposited config says in place of a machine's directories.
_PLACEHOLDER = "machine path removed by polyzymd study freeze: say where the runs are with data.local.yaml (polyzymd study locate)"


def _condition_configs(protocol: Any) -> list[str]:
    """Return the condition configs inside the study folder, relative to it."""
    return [
        str(path.relative_to(protocol.root))
        for path in protocol.conditions.values()
        if path.is_relative_to(protocol.root)
    ]


def without_machine_paths(text: str) -> str:
    """Return a config's text with the directories of one machine taken out, for the deposit.

    ``projects_directory`` and ``scratch_directory`` say where one machine
    keeps job files and runs; they become ``.`` and ``data``, with a comment
    saying how a reproducer points the study at their copy. The config hash
    leaves both out, so stored results still match. A ``Copied by polyzymd
    study init from <path>`` header keeps only the file name.
    """
    import re

    lines = []
    for line in text.splitlines():
        match = re.match(r"^(\s*)(projects_directory|scratch_directory):\s*.*$", line)
        if match:
            value = "." if match.group(2) == "projects_directory" else "data"
            line = f"{match.group(1)}{match.group(2)}: {value}  # {_PLACEHOLDER}"
        header = re.match(r"^# Copied by polyzymd study init from (.+)$", line)
        if header:
            line = f"# Copied by polyzymd study init from {Path(header.group(1).strip()).name}"
        lines.append(line)
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def _production_length_warnings(conditions: dict[str, Any]) -> list[str]:
    """Warn when the conditions' replicates were simulated for very different lengths.

    Uses the production length of each replicate in the manifest and the
    tolerance of :func:`polyzymd.analyses.study.production_length_warnings`.
    """
    from polyzymd.analyses.study import PRODUCTION_LENGTH_TOLERANCE

    lengths = {
        label: [
            r["production_ns"] for r in c.get("replicates", {}).values() if "production_ns" in r
        ]
        for label, c in conditions.items()
    }
    lengths = {label: v for label, v in lengths.items() if v}
    if len(lengths) < 2:
        return []
    longest = max(max(v) for v in lengths.values())
    shortest = min(min(v) for v in lengths.values())
    if longest <= 0 or (longest - shortest) / longest <= PRODUCTION_LENGTH_TOLERANCE:
        return []
    described = "; ".join(f"{label} {min(v):.4g}-{max(v):.4g} ns" for label, v in lengths.items())
    return [
        f"the conditions' production lengths differ ({described}); results compared across "
        f"them may reflect simulated time, so analyse with until {shortest:.4g}ns or extend the "
        "short runs, and say which in the methods"
    ]


def composition_warnings(label: str, config: Any, universe: Any) -> list[str]:
    """Return how the topology of ``label`` disagrees with what its config says was simulated.

    The residues that are neither protein, water nor ions are compared with
    the config: its substrate's ``residue_name`` should be among them, and
    the others are taken as polymer, which the config should enable. A
    disagreement means the deposited config does not describe the simulated
    system, which a reproducer would then build wrongly.
    """
    found = Counter(
        str(r).upper()
        for r in universe.select_atoms("not protein and not water").residues.resnames
        if str(r).upper() not in _IONS
    )
    notes: list[str] = []
    substrate = getattr(config, "substrate", None)
    substrate_name = str(substrate.residue_name).upper() if substrate is not None else None
    if substrate_name and substrate_name not in found:
        notes.append(
            f"{label}: the config names substrate residue {substrate_name}, which the topology "
            "does not contain"
        )
    others = {name: n for name, n in found.items() if name != substrate_name}
    polymers = getattr(config, "polymers", None)
    expects_polymer = bool(polymers is not None and getattr(polymers, "enabled", False))
    described = ", ".join(f"{name} {n}" for name, n in sorted(others.items()))
    if expects_polymer and not others:
        notes.append(
            f"{label}: the config enables polymers, but the topology has no residue besides "
            "protein, water, ions and the substrate"
        )
    elif not expects_polymer and others:
        notes.append(
            f"{label}: the topology contains residues {described} besides protein, water and "
            f"ions{' and the substrate' if substrate_name else ''}, but the config "
            f"{'has no substrate and ' if not substrate_name else ''}enables no polymers; the "
            "deposited config may not describe the simulated system"
        )
    return notes


def simulated_with(working_dir: Path) -> dict[str, Any]:
    """Return the software versions that built and ran a replicate, as the run recorded them.

    ``build`` comes from the replicate's ``build_manifest.json``, and
    ``segments`` lists each distinct (PolyzyMD, OpenMM, pixi environment)
    combination that ``progress.json`` records for its production segments,
    so a reproducer knows which engine produced the trajectories, not only
    which PolyzyMD analysed them.
    """
    from polyzymd.simulation.progress import load_progress

    found: dict[str, Any] = {}
    build = working_dir / "build_manifest.json"
    try:
        manifest = json.loads(build.read_text())
        found["build"] = {k: manifest.get(k) for k in ("polyzymd_version", "openmm_version")}
    except (OSError, ValueError):
        pass
    try:
        progress = load_progress(working_dir)
    except Exception:  # noqa: BLE001 - an unreadable progress file records nothing
        progress = None
    if progress is not None:
        combos = {
            (s.polyzymd_version, s.openmm_version, s.pixi_environment) for s in progress.segments
        }
        found["segments"] = [
            {"polyzymd_version": a, "openmm_version": b, "pixi_environment": c}
            for a, b, c in sorted(combos, key=lambda t: tuple(str(x) for x in t))
        ]
    return found


def _portable(value: Any, root: Path, key: str | None = None) -> Any:
    """Return ``value`` with every config path made relative to ``root``, or reduced to its name.

    A path outside the study folder is a fact about one machine, so only its
    file name is kept.
    """
    from polyzymd.config.loader import PATH_KEYS

    if isinstance(value, dict):
        return {k: _portable(v, root, k) for k, v in value.items()}
    if isinstance(value, list):
        return [_portable(v, root, key) for v in value]
    if key in PATH_KEYS and isinstance(value, str) and Path(value).is_absolute():
        path = Path(value)
        return str(path.relative_to(root)) if path.is_relative_to(root) else path.name
    return value


def _summary_row(label: str, index: int, universe: Any) -> dict[str, Any]:
    import numpy as np

    universe.trajectory[0]
    box = universe.dimensions if universe.dimensions is not None else [np.nan] * 6
    water = universe.select_atoms("water").residues
    residues = universe.residues
    ions = Counter(str(r) for r in residues.resnames if str(r).upper() in _IONS)
    other = Counter(
        str(r)
        for r in universe.select_atoms("not water and not protein").residues.resnames
        if str(r).upper() not in _IONS
    )
    return {
        "condition": label,
        "replicate": index,
        "atoms": len(universe.atoms),
        "box_a_A": round(float(box[0]), 3),
        "box_b_A": round(float(box[1]), 3),
        "box_c_A": round(float(box[2]), 3),
        "waters": len(water),
        "ions": "; ".join(f"{k} {v}" for k, v in sorted(ions.items())),
        "protein_residues": len(universe.select_atoms("protein").residues),
        "other_residues": "; ".join(f"{k} {v}" for k, v in sorted(other.items())),
    }


def _replicates(
    protocol: Any, deposit: Path, hashes: _Hashes, warnings: list[str]
) -> tuple[dict, list]:
    """Collect every replicate on this machine: hashes, engine inputs, final frame, summary row."""
    import warnings as python_warnings

    from polyzymd.analyses.study import Condition

    conditions: dict[str, Any] = {}
    rows: list[dict[str, Any]] = []
    for label, path in protocol.conditions.items():
        record: dict[str, Any] = {
            "config": str(path.relative_to(protocol.root))
            if path.is_relative_to(protocol.root)
            else str(path)
        }
        try:
            condition = Condition(
                label,
                path,
                protocol.equilibration,
                protocol.replicates,
                protocol.stride,
                protocol.data.get(label),
            )
        except ProtocolError as exc:
            warnings.append(
                f"{label}: the trajectories are not on this machine, so its replicates are not "
                f"hashed and its engine inputs and final frames are not deposited ({exc})"
            )
            conditions[label] = {**record, "replicates": {}}
            continue
        record["config_hash"] = condition.config_hash
        record["engine"] = str(getattr(condition.config, "engine", None) or "openmm")
        settings = condition.config.model_dump(mode="json")
        settings.get("output", {}).pop("projects_directory", None)
        settings.get("output", {}).pop("scratch_directory", None)
        record["resolved_config"] = _portable(settings, protocol.root)
        replicates: dict[str, Any] = {}
        for replicate in condition.replicates:
            with python_warnings.catch_warnings():
                python_warnings.simplefilter("ignore")
                replicate.universe()  # loading records the bond source in the provenance
            provenance = condition._provider.provenance_for(replicate.index)
            data_root = Path(condition.config.output.effective_scratch_directory)

            def relative(file: str | Path, data_root: Path = data_root) -> str:
                file = Path(file)
                return (
                    str(file.relative_to(data_root))
                    if file.is_relative_to(data_root)
                    else file.name
                )

            files = [
                {"path": relative(item.path), **hashes(Path(item.path))}
                for item in (provenance.topology, *provenance.trajectories)
            ]
            folder = (
                Path(DEPOSIT)
                / "engine_inputs"
                / condition_slug(label)
                / f"replicate_{replicate.index}"
            )
            inputs = []
            for source in _engine_inputs(provenance):
                target = _gzip_copy(source, protocol.root / folder / (source.name + ".gz"))
                inputs.append(
                    {
                        "path": str(target.relative_to(protocol.root)),
                        "source": source.name,
                        **hashes(target),
                    }
                )
            with python_warnings.catch_warnings():
                python_warnings.simplefilter("ignore")
                universe = replicate.universe()
                rows.append(_summary_row(label, replicate.index, universe))
                if replicate is condition.replicates[0]:
                    warnings.extend(composition_warnings(label, condition.config, universe))
                universe.trajectory[-1]
                final_pdb = (
                    protocol.root
                    / DEPOSIT
                    / "final_frames"
                    / condition_slug(label)
                    / f"replicate_{replicate.index}_final.pdb"
                )
                final_pdb.parent.mkdir(parents=True, exist_ok=True)
                universe.atoms.write(str(final_pdb))
            final = _gzip_copy(final_pdb, final_pdb.with_suffix(".pdb.gz"))
            final_pdb.unlink()
            dt = float(universe.trajectory.dt)
            replicates[str(replicate.index)] = {
                "files": files,
                "engine_inputs": inputs,
                "final_frame": {
                    "path": str(final.relative_to(protocol.root)),
                    "time_ps": float(universe.trajectory.time),
                    **hashes(final),
                },
                "production_frames": int(universe.trajectory.n_frames),
                "production_ns": round((universe.trajectory.n_frames - 1) * dt / 1000.0, 6),
                "frames_analysed": int(len(replicate.frames)),
                "first_analysed_ns": float(replicate.times[0]) if len(replicate.times) else None,
                "trajectory_variant": provenance.trajectory_variant,
                "bond_source": provenance.bond_source,
                "warnings": list(provenance.warnings),
                "simulated_with": simulated_with(Path(provenance.working_directory)),
            }
        conditions[label] = {**record, "replicates": replicates}
        for index, replicate_record in replicates.items():
            for text in replicate_record["warnings"]:
                warnings.append(f"{label} replicate {index}: {text}")
        unknown = [
            index
            for index, r in replicates.items()
            if not any(
                v.get("openmm_version")
                for v in [
                    r["simulated_with"].get("build", {}),
                    *r["simulated_with"].get("segments", []),
                ]
            )
        ]
        if unknown:
            warnings.append(
                f"{label}: replicates {', '.join(unknown)} record no OpenMM version (no "
                "build_manifest.json, and progress.json predates version recording); state the "
                "engine version in the methods"
            )
    return conditions, rows


def condition_slug(label: str) -> str:
    from polyzymd.analyses.study_scaffold import condition_folder

    return condition_folder(label)


def condition_restraints(protocol: Any) -> dict[str, list[dict[str, Any]]]:
    """Return each condition's enabled distance restraints, read from its config.

    These restraints are added to the system when it is built, so they act in
    every phase, production included: a restrained condition samples a
    biased ensemble. Equilibration-only position restraints are not listed.
    """
    from polyzymd.config.schema import SimulationConfig

    found: dict[str, list[dict[str, Any]]] = {}
    for label, path in protocol.conditions.items():
        try:
            config = SimulationConfig.from_yaml(path)
        except (OSError, ValueError):
            continue
        found[label] = [
            {
                "name": r.name,
                "type": r.type.value,
                "atom1": r.atom1.selection,
                "atom2": r.atom2.selection,
                "distance_A": r.distance,
                "force_constant_kJ_mol_nm2": r.force_constant,
            }
            for r in config.restraints
            if r.enabled
        ]
    return found


def _sampling(restraints: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """Answer checklist item 3c from the conditions' distance restraints."""
    restrained = {label: items for label, items in restraints.items() if items}
    if not restrained:
        return {"answer": "unbiased molecular dynamics: no condition has a distance restraint"}
    return {
        "answer": "restrained molecular dynamics: the listed distance restraints act in every "
        "phase, production included, so these conditions sample a biased ensemble; state the "
        "restraints and their purpose in the methods",
        "evidence": restrained,
    }


def _checklist(protocol: Any, manifest: dict[str, Any], meta: dict[str, Any]) -> dict[str, Any]:
    """Fill the Communications Biology checklist (2023) from the manifest; every answer is informational."""
    counts = {label: len(c["replicates"]) for label, c in manifest["conditions"].items()}
    restraints = condition_restraints(protocol)
    engines = sorted({c.get("engine") for c in manifest["conditions"].values() if c.get("engine")})
    resolved = [
        c["resolved_config"] for c in manifest["conditions"].values() if "resolved_config" in c
    ]
    first = resolved[0] if resolved else {}

    def item(answer: Any, evidence: Any = None) -> dict[str, Any]:
        return {"answer": answer, **({"evidence": evidence} if evidence is not None else {})}

    return {
        # An unsigned editorial, so it is cited by its title.
        "source": "Reliability and reproducibility checklist for molecular dynamics simulations "
        "(2023). Communications Biology 6:268. doi:10.1038/s42003-023-04653-0",
        "note": "Filled by polyzymd study freeze from manifest.json; informational. Review each answer.",
        "1a_equilibration_evidence": item(
            "per-replicate time series and detected equilibration starts are in each run's report",
            [f"results/{run}/report.json" for run in protocol.analyses],
        ),
        "1b_equilibration_and_production": item(
            f"equilibration window {protocol.equilibration} removed from every replicate; stride {protocol.stride}",
            {
                label: {r: v["frames_analysed"] for r, v in c["replicates"].items()}
                for label, c in manifest["conditions"].items()
            },
        ),
        "1c_replicates_and_statistics": item(
            "the replicate is the sampling unit; 95% Student t intervals and Welch tests with "
            "Benjamini-Hochberg correction",
            {
                "replicates_per_condition": counts,
                "at_least_3": all(n >= 3 for n in counts.values()),
            },
        ),
        "1d_independent_starting_configurations": item(
            "each replicate is built with its replicate number as the random seed", None
        ),
        "2a_connection_to_experiment": item(
            [e.get("description") or e.get("doi") for e in meta["related"]["experimental"]]
            or "TODO: describe"
        ),
        "3a_system_type": item(meta["system_type"] or "TODO: list (e.g. protein, polymer)"),
        "3b_model_accuracy": item(
            "TODO: justify the force field and water model", first.get("force_field")
        ),
        "3c_enhanced_sampling": _sampling(restraints),
        "4a_system_setup_table": item(SUMMARY),
        "4b_simulation_parameters": item(
            "thermodynamics, restraints, cutoffs, thermostat and barostat of each condition",
            {
                label: {
                    **{
                        k: c["resolved_config"].get(k)
                        for k in ("thermodynamics", "simulation_phases")
                    },
                    "restraints": restraints.get(label, []),
                }
                for label, c in manifest["conditions"].items()
                if "resolved_config" in c
            },
        ),
        "4c_software_versions": item(manifest["versions"], {"engines": engines}),
        "4d_coordinates_and_inputs": item(
            "initial structures in conditions/*/structures/, engine inputs and final frames in deposit/",
            None,
        ),
        "4e_custom_code_and_parameters": item(
            "analyses/ and figures/ hold the custom code; generated polymer parameters are in the "
            "serialized engine inputs in deposit/engine_inputs/"
        ),
    }


def _method(protocol: Any) -> str:
    import polyzymd

    runs = ", ".join(
        f"{run} ({entry.analysis or entry.function.qualname})"
        for run, entry in protocol.analyses.items()
    )
    return (
        f"Analysed with PolyzyMD {polyzymd.__version__}: {len(protocol.conditions)} conditions "
        f"({', '.join(protocol.conditions)}), equilibration window {protocol.equilibration}, "
        f"stride {protocol.stride}; analyses {runs or 'none'}. The replicate is the sampling unit."
    )


def _next_tag(root: Path) -> str:
    existing = (_git(root, "tag", "--list", "study-v*") or "").split()
    numbers = [int(t.split("v")[-1]) for t in existing if t.split("v")[-1].isdigit()]
    return f"study-v{max(numbers, default=0) + 1}"


def freeze(root: str | Path, *, tag: str | None = None) -> FreezeResult:
    """Freeze the study in ``root`` for publication; see the module docstring.

    Raises
    ------
    ProtocolError
        Only when the study file or its metadata cannot be read, or the tag
        already exists. Everything else is a warning in the result.
    """
    import yaml

    import polyzymd
    from polyzymd.analyses.study_file import load_study_file
    from polyzymd.analyses.study_git import git_state
    from polyzymd.analyses.study_metadata import (
        check_metadata,
        citation_cff,
        dump_cff,
        zenodo_json,
    )

    protocol = load_study_file(root)
    root = protocol.root
    meta, warnings = check_metadata(protocol.metadata)
    state = git_state(root)
    if state is None:
        warnings.append("the study is not a git repository, so freeze cannot commit or tag it")
    elif state["inputs_uncommitted"]:
        warnings.append(
            "uncommitted inputs are not part of the tagged study: "
            + ", ".join(state["inputs_uncommitted"])
        )
    tag = tag or (_next_tag(root) if state else None)
    if state and tag and _git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}"):
        raise ProtocolError(f"The tag {tag} already exists.", hint="Give another with --tag.")
    for run, why in stale_runs(protocol).items():
        warnings.append(f"run {run} may be stale: {'; '.join(why)}")

    deposit = root / DEPOSIT
    deposit.mkdir(exist_ok=True)
    hashes = _Hashes(deposit / _HASH_CACHE)
    conditions, rows = _replicates(protocol, deposit, hashes, warnings)
    hashes.save()
    for label, items in condition_restraints(protocol).items():
        conditions[label]["restraints"] = items
    warnings.extend(_production_length_warnings(conditions))

    released = date.today().isoformat()
    version = tag or "unversioned"
    if rows:
        with (root / SUMMARY).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    elif (root / SUMMARY).exists():
        warnings.append(
            f"{SUMMARY} was kept from an earlier freeze: no trajectories here to redo it"
        )
    else:
        warnings.append(f"no {SUMMARY}: no trajectories on this machine")

    from polyzymd.analyses.study_file import STUDY_FILE

    if state:
        tracked = (_git(root, "ls-files") or "").splitlines()
        results = (
            _git(root, "ls-files", "--others", "--exclude-standard", "--", "results") or ""
        ).splitlines()
        candidates = {p for p in (*tracked, *results) if p}
    else:
        candidates = {
            str(p.relative_to(root))
            for p in root.rglob("*")
            if p.is_file()
            and p.relative_to(root).parts[0] not in (DEPOSIT, ".git", "data.local.yaml")
        }
    study_files = sorted(p for p in candidates if p not in GENERATED and (root / p).is_file())
    manifest: dict[str, Any] = {
        "$schema": MANIFEST_SCHEMA_FILE,
        "schema": MANIFEST_SCHEMA,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tag": tag,
        "git": {
            "commit": state["commit"] if state else None,
            "inputs_uncommitted": state["inputs_uncommitted"] if state else None,
        },
        "study_file": {"path": STUDY_FILE, **hashes(protocol.path)},
        "versions": _versions(),
        "equilibration": protocol.equilibration,
        "stride": protocol.stride,
        "metadata": meta,
        "conditions": conditions,
        "analyses": {
            run: {
                "analysis": entry.analysis,
                "function": f"{entry.function.file.relative_to(root)}:{entry.function.qualname}"
                if entry.function
                else None,
                "settings": entry.settings or (entry.function.settings if entry.function else {}),
                "report": f"results/{run}/report.json",
            }
            for run, entry in protocol.analyses.items()
        },
        "trajectory_deposits": meta["related"]["trajectories"],
        "files": {p: hashes(root / p) for p in study_files},
        "cite": {
            "polyzymd": __import__("polyzymd.citation", fromlist=["citation_line"]).citation_line()
        },
        "deposit_without_machine_paths": [
            f"{config}: output.projects_directory and output.scratch_directory replaced, "
            "'Copied from' header reduced to a file name"
            for config in _condition_configs(protocol)
        ],
        "warnings": warnings,
    }
    hashes.save()
    (root / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    (root / CHECKLIST).write_text(
        "# The Communications Biology MD checklist, filled by polyzymd study freeze. Informational.\n"
        + yaml.safe_dump(_checklist(protocol, manifest, meta), sort_keys=False)
    )
    (root / CITATION).write_text(
        dump_cff(
            citation_cff(meta, version=version, released=released, commit=manifest["git"]["commit"])
        )
    )
    (root / ZENODO).write_text(
        json.dumps(
            zenodo_json(meta, version=version, released=released, method=_method(protocol)),
            indent=2,
        )
        + "\n"
    )
    gitignore = root / ".gitignore"
    lines = gitignore.read_text().splitlines() if gitignore.exists() else []
    added = False
    for entry, why in (
        (f"{DEPOSIT}/", "What polyzymd study freeze lays out for upload."),
        ("logs/", "Full logs of polyzymd commands; the console shows only warnings."),
    ):
        if entry not in lines:
            lines, added = [*lines, f"# {why}", entry], True
    if added:
        gitignore.write_text("\n".join([*lines, ""]))

    commit = None
    if state and tag:
        paths = [p for p in (*GENERATED, ".gitignore", "results") if (root / p).exists()]
        _git(root, "add", "--", *paths)
        if (
            _git(
                root,
                "commit",
                "--quiet",
                "-m",
                f"Freeze study as {tag} with polyzymd study freeze",
                "--",
                *paths,
            )
            is None
        ):
            warnings.append("git could not commit the frozen files (is user.name set?)")
        elif (
            _git(root, "tag", "-a", tag, "-m", f"Study frozen by PolyzyMD {polyzymd.__version__}")
            is None
        ):
            warnings.append(f"git could not create the tag {tag}")
        else:
            commit = (_git(root, "rev-parse", "HEAD") or "").strip() or None

    study_copy = deposit / "study"
    if study_copy.exists():
        shutil.rmtree(study_copy)
    study_copy.mkdir(parents=True)
    if commit:
        archive = subprocess.run(
            ["git", "-C", str(root), "archive", "--format=tar", tag], capture_output=True
        )
        subprocess.run(["tar", "-x", "-C", str(study_copy)], input=archive.stdout, check=False)
    else:
        for name in (*study_files, *(p for p in GENERATED if (root / p).exists())):
            target = study_copy / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / name, target)
    for config in _condition_configs(protocol):
        copied = study_copy / config
        if copied.is_file():
            copied.write_text(without_machine_paths(copied.read_text()))
    for name in (MANIFEST, CITATION, ZENODO):
        if (root / name).exists():
            shutil.copy2(root / name, deposit / name)
    shutil.copy2(
        Path(__file__).parent / "schemas" / MANIFEST_SCHEMA_FILE, deposit / MANIFEST_SCHEMA_FILE
    )
    from polyzymd.analyses.study_upload_guide import deposit_readme, prepare_upload

    # The deposit's README describes the study from its metadata; the study's
    # own README.md stays as written, inside the study.
    (deposit / "README.md").write_text(
        deposit_readme(
            study_name=root.name,
            tag=tag if commit else None,
            meta=meta,
            analyses=protocol.analyses,
            root=root,
        )
    )

    prepared = prepare_upload(
        deposit,
        study_name=root.name,
        tag=tag if commit else None,
        manifest=manifest,
        zenodo=json.loads((root / ZENODO).read_text()),
        warnings=warnings,
    )
    return FreezeResult(
        root,
        tag if commit else None,
        commit,
        deposit,
        manifest,
        warnings,
        guide=prepared["guide"],
        upload=prepared["upload"],
    )
