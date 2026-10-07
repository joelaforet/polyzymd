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
import json
import os
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
    #: True when git could not commit or tag; the deposit then names no tag.
    git_failed: bool = False


class _Hashes:
    """SHA-256 of files and their sizes, from the shared hash cache.

    The cache is :mod:`~polyzymd.analyses.shared.file_hashes`, so files that
    analyses already hashed are not read again.
    """

    def __call__(self, file: Path) -> dict[str, Any]:
        from polyzymd.analyses.shared.file_hashes import file_sha256

        return {"size": file.stat().st_size, "sha256": file_sha256(file)}


def _git(root: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), *arguments], capture_output=True, text=True, timeout=120
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def _versions(root: Path) -> dict[str, str | None]:
    """Return the versions of PolyzyMD, Python and the main packages, and the lock file's hash.

    OpenMM's is its full version, as replicate records give it. A package
    that reports version ``0.0.0`` (a conda build without its version in the
    metadata) has the version of its conda package record
    (``conda-meta/<name>-<version>-<build>.json`` of the environment), or
    ``None``. ``pixi.lock`` is the SHA-256 of ``environment/pixi.lock`` under
    ``root``, or else of the ``pixi.lock`` of the pixi workspace that runs
    freeze, which pins every package; ``None`` when neither exists.
    ``pixi.lock_file`` says which file it is, and whether it is deposited.
    """
    import hashlib
    import platform
    import sys

    import polyzymd
    from polyzymd.utils.version import get_openmm_version, pixi_workspace

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
            version = str(
                getattr(imported, "__version__", None) or getattr(imported, "version", None)
            )
            versions[module] = None if version == "0.0.0" else version
        except Exception:  # noqa: BLE001 - an absent or broken package is recorded as absent
            versions[module] = None
        if versions[module] is None:
            name = module.lower().replace(".", "-")
            for record in Path(sys.prefix, "conda-meta").glob(f"{name}-*.json"):
                try:
                    found = json.loads(record.read_text())
                except (OSError, ValueError):
                    continue
                if isinstance(found, dict) and found.get("name") == name:
                    versions[module] = found.get("version")
    versions["openmm"] = get_openmm_version()
    lock = root / "environment" / "pixi.lock"
    versions["pixi.lock_file"] = "environment/pixi.lock (deposited)"
    workspace = pixi_workspace()
    if not lock.is_file() and workspace:
        lock = workspace / "pixi.lock"
        versions["pixi.lock_file"] = (
            "pixi.lock of the pixi workspace that ran freeze (not deposited)"
        )
    if not lock.is_file():
        lock, versions["pixi.lock_file"] = None, None
    versions["pixi.lock"] = hashlib.sha256(lock.read_bytes()).hexdigest() if lock else None
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


def stale_runs(protocol: Any, conditions: dict[str, Any] | None = None) -> dict[str, list[str]]:
    """Return, for each analysis run, why its stored results may not match the study now.

    A run has one reason per difference between what produced its stored
    results and the study now: no stored results; a record's config hash,
    equilibration window, ``until`` (``until: common`` included, worked out
    from the runs when they are here), stride, function hash (with its
    folder's helper modules), or the content of a file it was given; the
    replicates found on disk; or the report's settings, selections, condition
    factors (which its trend tests used), comparison block or PolyzyMD
    version. Trajectories are read only to work out ``until: common``.
    ``conditions`` are the manifest's conditions (from
    :func:`_replicates`); with them, a replicate whose
    trajectory hashes differ from those its record names is a reason too.
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
    from polyzymd.analyses.study_file import entry_record, portable

    project_root = protocol.project.root if protocol.project is not None else None
    found_replicates = _replicates_on_disk(protocol)
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
        current_files = _named_files(
            entry.function.settings if entry.function is not None else entry.settings
        )
        common_end = _common_end(protocol, run)
        recorded_replicates: dict[str, set[int]] = {}
        for path in records:
            record = json.loads(path.read_text())
            label = record.get("condition")
            recorded_replicates.setdefault(label, set()).add(int(record.get("replicate", 0)))
            here = (conditions or {}).get(label, {}).get("replicates", {})
            on_disk = here.get(str(record.get("replicate")))
            if on_disk is not None:
                # The first file is the topology, the others the trajectories.
                now = sorted(item["sha256"] for item in on_disk["files"][1:])
                then = [item.get("sha256") for item in record.get("trajectories", [])]
                if None in then:
                    found.append(
                        f"no trajectory hash recorded for {label} replicate {record['replicate']}"
                    )
                elif now != sorted(then):
                    found.append(
                        f"the trajectories of {label} replicate {record['replicate']} changed"
                    )
            for name, sha in _recorded_files(record.get("arguments")).items():
                if name in current_files and current_files[name] != sha:
                    found.append(f"the content of {name} changed")
            if common_end is not None and record.get("until_ns") is not None:
                if abs(float(record["until_ns"]) - common_end) > 1e-6:
                    found.append(
                        f"until common is now {common_end:g}ns, not {record['until_ns']:g}ns"
                    )
            if label in hashes and hashes[label] and record.get("config_hash") != hashes[label]:
                found.append(f"the config of {label} changed")
            equilibration, until = protocol.window(run)
            if record.get("equilibration") != equilibration:
                found.append(f"equilibration {record.get('equilibration')} is not {equilibration}")
            if not _same_until(record.get("until_ns"), until):
                recorded = record.get("until_ns")
                found.append(
                    f"until {'none' if recorded is None else f'{recorded}ns'} is not "
                    f"{until or 'none'}"
                )
            if int(record.get("stride", 1)) != protocol.stride_of(run):
                found.append(f"stride {record.get('stride')} is not {protocol.stride_of(run)}")
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
        for label, replicates in found_replicates.items():
            if label in recorded_replicates and replicates != recorded_replicates[label]:
                found.append(
                    f"{label} has replicates {sorted(replicates)} on disk, the results "
                    f"{sorted(recorded_replicates[label])}"
                )
        study = provenance.get("study") or {}
        expected = entry.function.settings if entry.function is not None else entry.settings
        expected = portable(expected, protocol.root, project_root)
        ran_with = study.get("settings")
        if ran_with is not None:
            ran_with = portable(ran_with, protocol.root, project_root)
            for key in sorted(set(expected) | set(ran_with)):
                if _canon(expected.get(key)) != _canon(ran_with.get(key)):
                    found.append(f"setting {key} differs from the one the results were made with")
        if entry.function is not None and "selections" in study:
            if study["selections"] != dict(entry.function.selections):
                found.append("the selections differ from the ones the results were made with")
        if study and "entry" not in study:
            found.append(
                "the report predates the recording of its analysis entry; rerun polyzymd "
                "analyze (stored values are reused) so later changes to the entry are seen"
            )
        if study.get("entry") is not None:
            now_entry = entry_record(protocol, run)
            changed = sorted(
                key
                for key in set(now_entry) | set(study["entry"])
                if _canon(now_entry.get(key)) != _canon(study["entry"].get(key))
            )
            inner = []
            if "function" in changed and isinstance(now_entry.get("function"), dict):
                then = study["entry"].get("function") or {}
                inner = sorted(
                    key
                    for key in set(now_entry["function"]) | set(then)
                    if _canon(now_entry["function"].get(key)) != _canon(then.get(key))
                )
            for key in [k for k in changed if k != "function"] + inner:
                found.append(f"the entry's {key} changed since the report")
        reported = [c.get("label") for c in report.get("conditions", [])]
        listed = list(protocol.conditions)
        if reported and set(dict.fromkeys(reported)) != set(listed):
            missing_now = sorted(set(listed) - set(reported))
            extra = sorted(set(reported) - set(listed))
            found.append(
                "the report's conditions differ from study.yaml"
                + (f": not in the report {missing_now}" if missing_now else "")
                + (f"; no longer in study.yaml {extra}" if extra else "")
            )
        if "factors" in study:
            labels = [c.get("label") for c in report.get("conditions", [])]
            now = {k: v for k, v in protocol.factors.items() if k in labels}
            then = {k: v for k, v in study["factors"].items() if k in labels}
            if now != then:
                found.append("the condition factors changed since the report's trend tests")
            if study.get("comparison") != protocol.comparison:
                found.append("the comparison block changed since the report")
        if found:
            reasons[run] = sorted(set(found))
    return reasons


#: Pathspecs ``git add`` and ``git commit`` leave out when freezing: job files of
#: ``polyzymd analyze --submit``, logs, which name one machine's paths, and runs.
EXCLUDE_MACHINE_FILES = [
    ":(exclude)**/slurm/**",
    ":(exclude)**/slurm_logs/**",
    ":(exclude)**/logs/**",
    ":(exclude)**/runs/**",
]


#: Folders that hold one machine's job files, logs, runs or environments, never published.
MACHINE_FOLDERS = ("slurm", "slurm_logs", "logs", "runs")


def is_machine_file(path: str) -> bool:
    """Return whether ``path`` is never published: in a job, log, run or hidden folder, or hidden.

    Job and log folders (``slurm/``, ``slurm_logs/``, ``logs/``) name one
    machine's paths; ``runs/`` holds the simulations, which are deposited
    apart from the study; hidden ones (``.git``, ``.pixi``, ``.venv``) hold
    repositories and environments.
    """
    parts = Path(path).parts
    return any(part in MACHINE_FOLDERS for part in parts[:-1]) or any(
        part.startswith(".") and part not in (".gitignore", ".zenodo.json") for part in parts
    )


def drop_machine_files(folder: Path) -> None:
    """Delete every job, log, run and hidden folder under ``folder``, a copy being deposited."""
    for path in sorted(folder.rglob("*"), reverse=True):
        if path.is_dir() and (path.name in MACHINE_FOLDERS or path.name.startswith(".")):
            shutil.rmtree(path, ignore_errors=True)


def _canon(value: Any) -> str:
    """Return ``value`` as canonical JSON, with NaN and infinities as the strings reports use.

    A report stores ``missing: .nan`` as ``"NaN"``, and NaN never equals
    itself, so values are compared through this form.
    """
    import math

    def plain(item: Any) -> Any:
        if isinstance(item, float) and not math.isfinite(item):
            return "NaN" if math.isnan(item) else ("Infinity" if item > 0 else "-Infinity")
        if isinstance(item, dict):
            return {str(k): plain(v) for k, v in item.items()}
        if isinstance(item, (list, tuple)):
            return [plain(v) for v in item]
        return item

    return json.dumps(plain(value), sort_keys=True, default=str)


def committable(root: Path, paths: list[str]) -> list[str]:
    """Return the paths git can commit after ``git add``: the pathspecs that match a file.

    A folder holding only excluded files (a ``results/`` with only
    ``slurm/``) matches nothing, and naming it would make ``git commit`` fail.
    """
    return [
        p for p in paths if p.startswith(":(") or (_git(root, "ls-files", "--", p) or "").strip()
    ]


def report_problems(protocol: Any) -> list[str]:
    """Return one warning per stored report marked partial, naming its problems.

    A partial report covers only some conditions (``polyzymd analyze`` names
    the others in ``problems``); freezing it as if complete would hide that.
    """
    warnings = []
    for run in protocol.analyses:
        path = protocol.results_dir(run) / "report.json"
        try:
            report = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        if report.get("status", "complete") != "complete":
            warnings.append(
                f"run {run} has a {report.get('status')} report: "
                + "; ".join(report.get("problems") or ["no problems listed"])
            )
    return warnings


def _named_files(value: Any) -> dict[str, str]:
    """Return the SHA-256 of every existing file a setting names, by file name."""
    import hashlib

    files: dict[str, str] = {}
    if isinstance(value, dict):
        for item in value.values():
            files.update(_named_files(item))
    elif isinstance(value, list):
        for item in value:
            files.update(_named_files(item))
    elif isinstance(value, str) and len(value) < 4096 and Path(value).is_file():
        files[Path(value).name] = hashlib.sha256(Path(value).read_bytes()).hexdigest()
    return files


def _recorded_files(value: Any) -> dict[str, str]:
    """Return the SHA-256 of every file a record's arguments name, by file name."""
    files: dict[str, str] = {}
    if isinstance(value, dict):
        if set(value) == {"name", "sha256"}:
            return {value["name"]: value["sha256"]}
        for item in value.values():
            files.update(_recorded_files(item))
    elif isinstance(value, list):
        for item in value:
            files.update(_recorded_files(item))
    return files


def _replicates_on_disk(protocol: Any) -> dict[str, set[int]]:
    """Return each condition's replicates found on this machine, within ``replicates:``.

    A condition whose runs are not here is left out, so a study without its
    trajectories is not stale for it.
    """
    from polyzymd.analyses.study import with_data_dir
    from polyzymd.config.schema import SimulationConfig

    wanted = set(protocol.replicates) if protocol.replicates else None
    found: dict[str, set[int]] = {}
    for label, path in protocol.conditions.items():
        try:
            config = with_data_dir(SimulationConfig.from_yaml(path), protocol.data.get(label))
            indices = {int(i) for i, _ in config.discover_replicate_dirs()}
        except (OSError, ValueError):
            continue
        indices = indices & wanted if wanted is not None else indices
        if indices:
            found[label] = indices
    return found


def _common_end(protocol: Any, run: str) -> float | None:
    """Return the time ``until: common`` stands for now, in ns, or ``None``.

    ``None`` when the run's ``until`` is not ``common`` or the runs are not
    here to work it out.
    """
    if protocol.window(run)[1] != "common":
        return None
    from polyzymd.analyses.study import Study

    try:
        study = Study.from_configs(
            dict(protocol.conditions),
            equilibration=protocol.window(run)[0],
            replicates=protocol.replicates,
            stride=protocol.stride_of(run),
            data=dict(protocol.data),
            until="common",
        )
        return next(iter(study)).until_ns
    except Exception:  # noqa: BLE001 - without the runs, until: common cannot be checked
        return None


def _own_windows(protocol: Any) -> str:
    """Name the analyses that set their own window, for the checklist and methods text."""
    own = []
    for run, entry in protocol.analyses.items():
        if entry.equilibration is None and entry.until is None and entry.stride is None:
            continue
        equilibration, until = protocol.window(run)
        own.append(
            f"{run} {equilibration}"
            + (f" until {until}" if until else "")
            + (f" stride {entry.stride}" if entry.stride else "")
        )
    return f" (analyses with their own window: {', '.join(own)})" if own else ""


def _same_until(recorded_ns: Any, until: str | None) -> bool:
    """Return whether a record's ``until_ns`` is the window end ``until`` (None: no end)."""
    if until is None or recorded_ns is None:
        return until is None and recorded_ns is None
    if until == "common":
        # Compared with the end worked out from the runs, by _common_end.
        return True
    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    value, unit = parse_time_string(until)
    return abs(convert_time(value, unit, "ns") - float(recorded_ns)) < 1e-9


def _gzip_copy(source: Path, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as raw, gzip.open(target, "wb") as packed:
        shutil.copyfileobj(raw, packed)
    return target


def _engine_inputs(provenance: Any, config: Any = None) -> list[Path]:
    """Return the engine input files of a replicate: OpenMM system XML and topology, or GROMACS files.

    GROMACS files are those PolyzyMD writes, by name
    (:func:`~polyzymd.analyses.shared.gromacs.run_input_files`): another
    ``.top`` or ``.mdp`` left in the run folder is not deposited.
    """
    from polyzymd.analyses.shared.gromacs import run_input_files
    from polyzymd.analyses.shared.loader import openmm_system_file

    topology = Path(provenance.topology.path)
    files: list[Path] = []
    if (provenance.config_engine or "openmm") == "gromacs":
        files.extend(run_input_files(topology.parent, config))
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
_PLACEHOLDER = "machine path removed by polyzymd: say where the runs are with data.local.yaml (polyzymd study locate)"


def deposited_entry(file: Path, root: Path, configs: set[Path]) -> dict[str, Any]:
    """Return the size and SHA-256 of ``file`` as freeze deposits it from the folder ``root``.

    A condition config (one of ``configs``, resolved) is deposited without
    machine paths (:func:`without_machine_paths`), so its entry describes
    that text; any other file is deposited as it is.
    """
    import hashlib

    if file.resolve() in configs:
        text = file.read_bytes().decode()
        deposited = without_machine_paths(text, file.parent, root)
        if deposited != text:
            data = deposited.encode()
            return {"size": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    return _Hashes()(file)


def _condition_configs(protocol: Any) -> list[str]:
    """Return the condition configs inside the study folder, relative to it."""
    return [
        str(path.relative_to(protocol.root))
        for path in protocol.conditions.values()
        if path.is_relative_to(protocol.root)
    ]


def without_machine_paths(text: str, folder: Path | None = None, root: Path | None = None) -> str:
    """Return a config's text with the directories of one machine taken out, for the deposit.

    An absolute ``projects_directory`` or ``scratch_directory`` says where one
    machine keeps job files and runs; it becomes ``.`` or ``data``, with a
    comment saying how a reproducer points the study at their copy. Relative
    values stay. The config hash leaves both out, so stored results still
    match. A ``Copied by polyzymd study init from <path>`` header keeps only
    the file name. With ``folder`` (the folder that holds the config) and
    ``root`` (the study or project folder), an absolute input path
    (``config.loader.PATH_KEYS``) inside ``root`` becomes a path relative to
    ``folder``.

    The lines are rewritten in place, so comments stay. When an input path
    changes, either output key keeps an absolute value, or the text no
    longer reads as YAML (a flow-style ``output: {...}`` or a block scalar),
    the config is read and written back with the changes, without its
    comments. A config with no machine path is returned unchanged.
    """
    import yaml

    rewritten = _rewrite_machine_lines(text)
    data = (yaml.safe_load(text) or {}) if folder is not None and root is not None else {}
    moved = [
        (container, key, Path(value).resolve())
        for container, key, value in _input_paths(data)
        if Path(value).is_absolute() and Path(value).resolve().is_relative_to(root.resolve())
    ]
    for container, key, path in moved:
        container[key] = os.path.relpath(path, folder.resolve())
    if not moved:
        try:
            output = (yaml.safe_load(rewritten) or {}).get("output") or {}
            if not any(_is_machine_path(output.get(key)) for key in _MACHINE_KEYS):
                return rewritten
        except (yaml.YAMLError, AttributeError):
            pass
        data = yaml.safe_load(text) or {}
    output = data.get("output")
    if isinstance(output, dict):
        for key, value in _MACHINE_KEYS.items():
            if _is_machine_path(output.get(key)):
                output[key] = value
    header = [line for line in rewritten.splitlines()[:1] if line.startswith("# Copied by")]
    return "\n".join([*header, f"# {_PLACEHOLDER}", yaml.safe_dump(data, sort_keys=False)])


def _input_paths(value: Any, key: str | None = None) -> list[tuple[Any, Any, str]]:
    """Return ``(container, key, path)`` for every input path (``PATH_KEYS``) in a config's data."""
    from polyzymd.config.loader import PATH_KEYS

    pairs = (
        value.items()
        if isinstance(value, dict)
        else enumerate(value) if isinstance(value, list) else ()
    )
    found = []
    for index, item in list(pairs):
        name = index if isinstance(value, dict) else key
        if name in PATH_KEYS and isinstance(item, str):
            found.append((value, index, item))
        else:
            found += _input_paths(item, name)
    return found


def _outside_inputs(protocol: Any) -> list[str]:
    """Warn for each absolute input path of a condition config outside the study and its project."""
    import yaml

    inside = [protocol.root] + ([protocol.project.root] if protocol.project is not None else [])
    warnings = []
    for label, config in protocol.conditions.items():
        try:
            data = yaml.safe_load(config.read_text()) or {}
        except (OSError, yaml.YAMLError):
            continue
        for _, _, value in _input_paths(data):
            path = Path(value)
            if path.is_absolute() and not any(
                path.resolve().is_relative_to(folder.resolve()) for folder in inside
            ):
                # Only the file name: the warning is published, the path names a machine.
                warnings.append(
                    f"{label}: the config names {path.name} outside the study, so the deposit "
                    "does not hold it; copy it into the study with polyzymd study add-condition"
                )
    return warnings


#: What an absolute output directory becomes in a deposited config.
_MACHINE_KEYS = {"projects_directory": ".", "scratch_directory": "data"}


def _is_machine_path(value: Any) -> bool:
    """Return whether a config value is an absolute path, starts with ``~`` or names a ``$VAR``."""
    return isinstance(value, str) and (
        value.startswith("~") or "$" in value or Path(value).is_absolute()
    )


def _rewrite_machine_lines(text: str) -> str:
    """Rewrite the absolute ``projects_directory``/``scratch_directory`` lines and the copy header."""
    import re

    import yaml

    lines = []
    for full in text.splitlines(keepends=True):
        line = full.rstrip("\r\n")
        end = full[len(line) :]
        match = re.match(r"^(\s*)(projects_directory|scratch_directory):(.*)$", line)
        try:
            value = yaml.safe_load(match.group(3)) if match else None
        except yaml.YAMLError:
            value = match.group(3).strip()
        if _is_machine_path(value):
            line = f"{match.group(1)}{match.group(2)}: {_MACHINE_KEYS[match.group(2)]}  # {_PLACEHOLDER}"
        header = re.match(r"^# Copied by polyzymd(?: study init)? from (.+)$", line)
        if header:
            line = f"# Copied by polyzymd from {Path(header.group(1).strip()).name}"
        lines.append(line + end)
    return "".join(lines)


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
    # Co-solvents (a surfactant, DMSO) are named by their residue_name, or
    # the first three letters of their name, as the builder names them.
    cosolvents = {
        str(cs.residue_name or cs.name[:3]).upper()
        for cs in getattr(getattr(config, "solvent", None), "co_solvents", None) or []
    }
    others = {
        name: n for name, n in found.items() if name != substrate_name and name not in cosolvents
    }
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
            f"{label}: the topology contains residues {described} besides protein, water, "
            f"ions, co-solvents{' and the substrate' if substrate_name else ''}, but the config "
            f"{'has no substrate and ' if not substrate_name else ''}enables no polymers; the "
            "deposited config may not describe the simulated system"
        )
    return notes


def _hashes_recorded(condition: Any, replicate: int, provenance: Any) -> bool:
    """Return whether the run recorded the hash of every trajectory file the replicate reads.

    Read through the condition's simulation engine, from the hashes the runner
    recorded and those ``polyzymd hash-trajectories`` recorded. Trajectories
    loaded without a simulation engine count as recorded.
    """
    recorded = condition._provider.recorded_trajectory_hashes(replicate)
    if recorded is None:
        return True
    paths = [Path(item.path).resolve() for item in provenance.trajectories]
    return bool(paths) and all(path in recorded for path in paths)


def _missing_build_files(working_dir: Path) -> list[str]:
    """Return the files ``build_manifest.json`` lists that are not in the run directory."""
    try:
        manifest = json.loads((working_dir / "build_manifest.json").read_text())
    except (OSError, ValueError):
        return []
    return sorted(
        name for name in (manifest.get("artifacts") or {}) if not (working_dir / name).exists()
    )


def simulated_with(working_dir: Path) -> dict[str, Any]:
    """Return the software versions that built and ran a replicate, as the run recorded them.

    ``build`` comes from the replicate's ``build_manifest.json``, and
    ``segments`` lists each distinct (PolyzyMD, OpenMM, pixi environment,
    OpenMM platform and properties) combination that ``progress.json``
    records for its production segments,
    so a reproducer knows which engine produced the trajectories, not only
    which PolyzyMD analysed them. For a GROMACS run, ``gromacs_version`` is
    the version that ``gmx mdrun`` wrote into ``gromacs/prod.log``.
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
            (
                s.polyzymd_version,
                s.openmm_version,
                s.pixi_environment,
                json.dumps(s.openmm_platform, sort_keys=True),
            )
            for s in progress.segments
        }
        found["segments"] = [
            {
                "polyzymd_version": a,
                "openmm_version": b,
                "pixi_environment": c,
                "openmm_platform": json.loads(d),
            }
            for a, b, c, d in sorted(combos, key=lambda t: tuple(str(x) for x in t))
        ]
    from polyzymd.engines.gromacs.engine import GromacsEngine

    try:
        with open(working_dir / GromacsEngine.engine_subdir / "prod.log", errors="ignore") as log:
            versions = {
                line.split(":", 1)[1].strip() for line in log if line.startswith("GROMACS version:")
            }
    except OSError:
        versions = set()
    if versions:
        found["gromacs_version"] = ", ".join(sorted(versions))
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
    # Engine inputs and final frames are recorded by their path in the deposit
    # and its zips, where a project puts each study's in a folder of its own.
    in_deposit = Path(protocol.project_label or "")
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
        except ProtocolError:
            # The error names this machine's paths, which the manifest must not.
            warnings.append(
                f"{label}: the trajectories are not on this machine, so its replicates are not "
                "hashed and its engine inputs and final frames are not deposited; run polyzymd "
                "study check on a machine that has them"
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
            for source in _engine_inputs(provenance, condition.config):
                target = _gzip_copy(source, protocol.root / folder / (source.name + ".gz"))
                inputs.append(
                    {
                        "path": str(
                            "engine_inputs" / in_deposit / target.relative_to(deposit / "engine_inputs")
                        ),
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
            replicates[str(replicate.index)] = {
                "files": files,
                "engine_inputs": inputs,
                "final_frame": {
                    "path": str(
                        "final_frames" / in_deposit / final.relative_to(deposit / "final_frames")
                    ),
                    "time_ps": float(universe.trajectory.time),
                    **hashes(final),
                },
                "production_frames": int(universe.trajectory.n_frames),
                "production_ns": round(replicate.production_ns, 6),
                "trajectory_variant": provenance.trajectory_variant,
                "bond_source": provenance.bond_source,
                "warnings": list(provenance.warnings),
                "simulated_with": simulated_with(Path(provenance.working_directory)),
                "missing_build_files": _missing_build_files(Path(provenance.working_directory)),
                "hashes_recorded": _hashes_recorded(condition, replicate.index, provenance),
            }
        conditions[label] = {**record, "replicates": replicates}
        for index, replicate_record in replicates.items():
            for text in replicate_record["warnings"]:
                warnings.append(f"{label} replicate {index}: {text}")
        engine = record.get("engine", "openmm")
        unknown = [
            index
            for index, r in replicates.items()
            if engine == "openmm"
            and not any(
                v.get("openmm_version")
                for v in [
                    r["simulated_with"].get("build", {}),
                    *r["simulated_with"].get("segments", []),
                ]
            )
        ]
        unhashed = [index for index, r in replicates.items() if not r["hashes_recorded"]]
        if unhashed:
            # From a project's folder the study is named; from the study's, it is ".".
            where = (
                f"--study {protocol.root.relative_to(protocol.project.root)} in the project folder"
                if protocol.project is not None
                else "--study . in the study folder"
            )
            warnings.append(
                f"{label}: replicates {', '.join(unhashed)} have no recorded trajectory "
                f"hashes; record them, once, by running polyzymd hash-trajectories {where}, "
                "so anyone can check the trajectories without hashing them again (it writes "
                "trajectory_hashes.json beside the trajectories, so whoever can write there "
                "runs it)"
            )
        if unknown:
            warnings.append(
                f"{label}: replicates {', '.join(unknown)} record no OpenMM version (no "
                "build_manifest.json, and progress.json predates version recording); state the "
                "engine version in the methods"
            )
        no_gromacs = [
            index
            for index, r in replicates.items()
            if engine == "gromacs" and not r["simulated_with"].get("gromacs_version")
        ]
        if no_gromacs:
            warnings.append(
                f"{label}: replicates {', '.join(no_gromacs)} record no GROMACS version (no "
                "gromacs/prod.log in the run directory); state the engine version in the methods"
            )
        for index, replicate_record in replicates.items():
            for name in replicate_record.get("missing_build_files", []):
                warnings.append(
                    f"{label} replicate {index}: build_manifest.json lists {name}, which is "
                    "missing from the run directory"
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
    polymers = _has_polymers(protocol)

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
            f"equilibration window {protocol.equilibration} removed from every replicate"
            + _own_windows(protocol)
            + f"; stride {protocol.stride}",
            # Each analysis has its own window, so its report holds the frames it used.
            {run: f"results/{run}/report.json: frames_per_replicate" for run in protocol.analyses},
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
            "each replicate's starting structure is built with its replicate number as the "
            "seed of Packmol" + (" and of polymer draws" if polymers else "") + "; the replicate "
            "number also seeds the initial velocities and the thermostat noise of each stage",
            None,
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
                    # The config names OpenMM's algorithms; GROMACS runs these instead.
                    **(
                        {"gromacs_production": _gromacs_production(c["resolved_config"])}
                        if c.get("engine") == "gromacs"
                        else {}
                    ),
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
            "analyses/ and figures/ hold the custom code; the parameters of every molecule"
            + (", generated polymers included," if polymers else "")
            + " are in the serialized engine inputs in deposit/engine_inputs/"
        ),
    }


def _gromacs_production(config: dict[str, Any]) -> dict[str, str]:
    """Return the integrator and coupling settings GROMACS runs for the production phase."""
    from polyzymd.exporters.gromacs import BAROSTAT_MAP, THERMOSTAT_MAP

    production = config["simulation_phases"]["production"]
    integrator, tcoupl = THERMOSTAT_MAP[production["thermostat"]]
    barostat = production.get("barostat") if production["ensemble"] == "NPT" else None
    pcoupl, pcoupltype = BAROSTAT_MAP.get(barostat, ("no", "isotropic"))
    return {"integrator": integrator, "tcoupl": tcoupl, "pcoupl": pcoupl, "pcoupltype": pcoupltype}


def _has_polymers(protocol: Any) -> bool:
    """Return whether any condition's config enables polymers."""
    from polyzymd.config.schema import SimulationConfig

    for path in protocol.conditions.values():
        try:
            polymers = SimulationConfig.from_yaml(path).polymers
        except Exception:  # noqa: BLE001 - an unreadable config is reported elsewhere
            continue
        if polymers is not None and getattr(polymers, "enabled", False):
            return True
    return False


def _method(protocol: Any) -> str:
    import polyzymd

    runs = ", ".join(
        f"{run} ({entry.analysis or entry.function.qualname})"
        for run, entry in protocol.analyses.items()
    )
    return (
        f"Analysed with PolyzyMD {polyzymd.__version__}: {len(protocol.conditions)} conditions "
        f"({', '.join(protocol.conditions)}), equilibration window {protocol.equilibration}"
        f"{_own_windows(protocol)}, stride {protocol.stride}; analyses {runs or 'none'}. The "
        "replicate is the sampling unit."
    )


def _next_tag(root: Path, prefix: str = "study") -> str:
    """Return ``<prefix>-v<n>``, with ``n`` one more than the highest existing such tag."""
    existing = (_git(root, "tag", "--list", f"{prefix}-v*") or "").split()
    numbers = [int(t.split("v")[-1]) for t in existing if t.split("v")[-1].isdigit()]
    return f"{prefix}-v{max(numbers, default=0) + 1}"


def group_warnings(warnings: list[str], names: list[str]) -> list[str]:
    """Return ``warnings`` with those that differ only in a leading name merged into one.

    A warning ``"<name>: <text>"`` for several names in ``names`` (conditions,
    or ``"<study>: <condition>"`` in a project) becomes ``"<name>, <name>:
    <text>"``, at the place of its first one. Nothing is dropped; the same
    text for eight conditions is one line instead of eight.
    """
    order: list[str] = []
    grouped: dict[str, list[str]] = {}
    for warning in warnings:
        name = next(
            (n for n in sorted(names, key=len, reverse=True) if warning.startswith(f"{n}: ")),
            None,
        )
        if name is None:
            order.append(warning)
            continue
        text = warning[len(name) + 2 :]
        if text not in grouped:
            grouped[text] = []
            order.append("\0" + text)
        grouped[text].append(name)
    return [
        f"{', '.join(grouped[item[1:]])}: {item[1:]}" if item.startswith("\0") else item
        for item in order
    ]


def _has_identity(root: Path) -> bool:
    """Return whether git has a user name and email for commits in ``root``, set or in the environment."""
    return all(
        os.environ.get(f"GIT_COMMITTER_{key.upper()}")
        or (_git(root, "config", f"user.{key}") or "").strip()
        for key in ("name", "email")
    )


def _git_preflight(
    root: Path, tag: str | None, what: str
) -> tuple[dict[str, Any] | None, str | None, list[str]]:
    """Return the git state of ``root``, the tag to freeze it as, and warnings about both.

    ``what`` is ``"study"`` or ``"project"``: it names the folder in the
    warnings and prefixes the default tag, the next ``<what>-v<n>``. The
    warnings say when ``root`` is not a git repository and which inputs are
    uncommitted.

    Raises
    ------
    ProtocolError
        If the tag already exists.
    """
    from polyzymd.analyses.study_git import git_state

    state = git_state(root)
    warnings = []
    if state is None:
        warnings.append(f"the {what} is not a git repository, so freeze cannot commit or tag it")
    elif state["inputs_uncommitted"]:
        # The tag and the deposit hold the committed files; the manifest must
        # describe exactly those, so every input is committed first.
        listed = state["inputs_uncommitted"]
        raise ProtocolError(
            f"The {what} has {len(listed)} uncommitted input files, which the tag and the "
            f"deposit would not contain: {', '.join(listed[:10])}"
            + (f" and {len(listed) - 10} more" if len(listed) > 10 else "")
            + ".",
            hint=f"Commit them first: git -C {root} add -A && git -C {root} commit -m 'Inputs "
            "for publication'. List files that must stay private in .gitignore. Results, logs "
            "and the deposit are committed by freeze itself.",
        )
    elif not _has_identity(root):
        raise ProtocolError(
            "git has no user name and email here, so freeze could not commit and tag.",
            hint=f"Set them: git -C {root} config user.name 'Your Name' && git -C {root} config "
            "user.email you@example.org",
        )
    tag = tag or (_next_tag(root, what) if state else None)
    if state and tag and _git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}"):
        raise ProtocolError(f"The tag {tag} already exists.", hint="Give another with --tag.")
    return state, tag, warnings


def deposited_folders() -> tuple[str, ...]:
    """Return the folders whose files freeze deposits: those ``study init`` and ``project init`` make."""
    from polyzymd.analyses.project_scaffold import PROJECT_FOLDERS
    from polyzymd.analyses.study_scaffold import FOLDERS

    return tuple(dict.fromkeys([*FOLDERS, *PROJECT_FOLDERS]))

#: Files at the top of a study or project that freeze deposits, besides those it writes.
DEPOSITED_FILES = ("study.yaml", "project.yaml", "data.example.yaml", ".gitignore")
#: Prefixes of other top-level files that freeze deposits (README.md, LICENSE, ...).
DEPOSITED_PREFIXES = ("README", "LICENSE")


def is_deposited_name(root: Path, path: str) -> bool:
    """Return whether freeze deposits ``path``, relative to the study or project ``root``.

    Freeze deposits only names that PolyzyMD chooses: ``study.yaml``,
    ``project.yaml``, ``data.example.yaml``, ``README*``, ``LICENSE*``,
    ``.gitignore``, the files freeze writes, and the files under the folders
    ``study init`` and ``project init`` make (:func:`deposited_folders`: ``conditions/``,
    ``structures/``, ``analyses/``, ``figures/``, ``results/``,
    ``environment/``, ``stats/`` ...). In a project, the same
    rule applies inside each study folder (a folder that holds ``study.yaml``).
    """
    parts = Path(path).parts
    if len(parts) > 1 and (root / parts[0] / "study.yaml").is_file():
        parts = parts[1:]
    if len(parts) > 1:
        return parts[0] in deposited_folders()
    name = parts[0]
    return (
        name in DEPOSITED_FILES
        or name in GENERATED
        or name.startswith(DEPOSITED_PREFIXES)
    )


#: Untracked paths that freeze lists in a git repository; it commits them.
_UNTRACKED = ("results", ".gitignore")


def _candidate_files(root: Path, state: dict[str, Any] | None) -> list[str]:
    """Return the files under ``root`` that freeze may publish, before the name rule.

    In a git repository (``state`` given), the tracked files, the untracked
    files under ``results/`` and an untracked ``.gitignore`` (freeze writes
    one and commits it); otherwise every file. Job, log and hidden files
    (:func:`is_machine_file`), anything under a ``deposit/`` folder and
    ``data.local.yaml`` are left out.
    """
    if state is not None:
        # -z: names separated by NUL and not quoted, as git quotes non-ASCII names.
        listed = (_git(root, "ls-files", "-z") or "").split("\0") + (
            _git(root, "ls-files", "-z", "--others", "--exclude-standard", "--", *_UNTRACKED)
            or ""
        ).split("\0")
    else:
        listed = [str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()]
    return sorted(
        {
            p
            for p in listed
            if p
            and DEPOSIT not in Path(p).parts
            and Path(p).name != "data.local.yaml"
            and not is_machine_file(p)
            and (root / p).is_file()
        }
    )


def _listed_files(root: Path, state: dict[str, Any] | None) -> list[str]:
    """Return the files under ``root`` that freeze hashes and deposits, relative to it.

    The files of :func:`_candidate_files` whose names :func:`is_deposited_name`
    accepts. Other files, such as notes or a copied trajectory, stay out of
    the deposit; :func:`left_out_files` names them.
    """
    return [p for p in _candidate_files(root, state) if is_deposited_name(root, p)]


def left_out_files(root: Path, state: dict[str, Any] | None) -> str | None:
    """Return a warning that names the files freeze does not deposit, or ``None``.

    A folder is named once, with a trailing ``/``, for all the files in it.
    """
    entries = []
    for path in _candidate_files(root, state):
        if is_deposited_name(root, path):
            continue
        parts = Path(path).parts
        study = len(parts) > 1 and (root / parts[0] / "study.yaml").is_file()
        depth = 2 if study else 1
        entry = "/".join(parts[:depth]) + ("/" if len(parts) > depth else "")
        if entry not in entries:
            entries.append(entry)
    if not entries:
        return None
    return (
        f"not deposited: {', '.join(entries)}. Freeze deposits only study.yaml, project.yaml, "
        "data.example.yaml, README*, LICENSE*, the files it writes, and "
        f"{', '.join(f'{f}/' for f in deposited_folders())}"
        "; move a file there to publish it"
    )


#: Entries freeze adds to ``.gitignore``, each with the comment above it.
_IGNORED = (
    (f"{DEPOSIT}/", "What polyzymd freeze lays out for upload."),
    ("logs/", "Full logs of polyzymd commands; the console shows only warnings."),
    ("data.local.yaml", "Where this machine keeps the trajectories."),
    ("runs/", "Simulation runs: trajectories never go into git or the deposit."),
    (".polymer_cache/", "Fragments and chains that dynamic polymer builds write."),
)


def _write_citation(
    root: Path,
    meta: dict[str, Any],
    *,
    version: str,
    released: str,
    commit: str | None,
    method: str,
) -> None:
    """Write ``CITATION.cff`` and ``.zenodo.json`` from ``meta``."""
    from polyzymd.analyses.study_metadata import citation_cff, dump_cff, zenodo_json

    (root / CITATION).write_text(
        dump_cff(citation_cff(meta, version=version, released=released, commit=commit))
    )
    (root / ZENODO).write_text(
        json.dumps(zenodo_json(meta, version=version, released=released, method=method), indent=2)
        + "\n"
    )


def _write_gitignore(root: Path) -> None:
    """Add the entries of :data:`_IGNORED` that ``root/.gitignore`` lacks, creating it if needed."""
    gitignore = root / ".gitignore"
    lines = gitignore.read_text().splitlines() if gitignore.exists() else []
    missing = [(entry, why) for entry, why in _IGNORED if entry not in lines]
    if missing:
        lines += [line for entry, why in missing for line in (f"# {why}", entry)]
        gitignore.write_text("\n".join([*lines, ""]))


def _commit_and_tag(
    root: Path, paths: list[str], tag: str, what: str, warnings: list[str]
) -> str | None:
    """Commit ``paths`` of ``root`` and tag the commit; return it, or ``None`` with a warning.

    Job and log files under the paths are never committed
    (:data:`EXCLUDE_MACHINE_FILES`). ``what`` (``"study"`` or ``"project"``)
    names the folder in the commit and tag messages.
    """
    import polyzymd

    paths = [*paths, *EXCLUDE_MACHINE_FILES]
    _git(root, "add", "--", *paths)
    paths = committable(root, paths)
    message = f"Freeze {what} as {tag} with polyzymd {what} freeze"
    if _git(root, "commit", "--quiet", "-m", message, "--", *paths) is None:
        warnings.append("git could not commit the frozen files (is user.name set?)")
        return None
    note = f"{what.capitalize()} frozen by PolyzyMD {polyzymd.__version__}"
    if _git(root, "tag", "-a", tag, "-m", note) is None:
        warnings.append(f"git could not create the tag {tag}")
        return None
    return (_git(root, "rev-parse", "HEAD") or "").strip() or None


def _drop_tag(
    root: Path, manifest: dict[str, Any], meta: dict[str, Any], released: str, method: str
) -> None:
    """Write the manifest and citation files again without a tag, after git could not make it."""
    manifest["tag"] = None
    _write_citation(
        root, meta, version="unversioned", released=released, commit=None, method=method
    )
    for name in (CITATION, ZENODO):
        if name in manifest["files"]:
            manifest["files"][name] = _Hashes()(root / name)
    (root / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")


def _copy_frozen_folder(
    root: Path,
    deposit: Path,
    tag: str | None,
    commit: str | None,
    files: list[str],
    configs: list[str],
) -> None:
    """Write ``deposit/study``: the files ``files`` of the frozen folder.

    They come from the tagged commit (``git archive``) or, without a commit,
    from the folder. The condition configs ``configs`` (relative to
    ``root``) are written there without machine paths
    (:func:`without_machine_paths`).
    """
    copy = deposit / "study"
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    if commit:
        tracked = set((_git(root, "ls-tree", "-r", "-z", "--name-only", tag) or "").split("\0"))
        names = [name for name in files if name in tracked]
        # Literal pathspecs, so a[1].csv is a name and not a pattern; in chunks,
        # to keep each command line short. Without paths, git archive would
        # write every tracked file.
        for start in range(0, len(names), 500):
            archive = subprocess.run(
                [
                    *("git", "-C", str(root), "--literal-pathspecs", "archive", "--format=tar"),
                    *(tag, "--", *names[start : start + 500]),
                ],
                capture_output=True,
            )
            subprocess.run(["tar", "-x", "-C", str(copy)], input=archive.stdout, check=False)
    else:
        for name in files:
            target = copy / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / name, target)
    drop_machine_files(copy)
    for config in configs:
        copied = copy / config
        if copied.is_file():
            # Bytes, so that CRLF line ends stay as the manifest hashed them.
            text = copied.read_bytes().decode()
            deposited = without_machine_paths(text, (root / config).parent, root)
            if deposited != text:
                copied.write_bytes(deposited.encode())


def _finish_deposit(
    root: Path,
    deposit: Path,
    tag: str | None,
    commit: str | None,
    manifest: dict[str, Any],
    readme: str,
    warnings: list[str],
    kind: str = "study",
) -> FreezeResult:
    """Copy the manifest, citation files and manifest schema into ``deposit``, write its README and upload guide.

    The upload folder and ``UPLOAD.md`` come from
    :func:`~polyzymd.analyses.study_upload_guide.prepare_upload`.
    """
    from polyzymd.analyses.study_upload_guide import prepare_upload

    for name in (MANIFEST, CITATION, ZENODO):
        if (root / name).exists():
            shutil.copy2(root / name, deposit / name)
    shutil.copy2(
        Path(__file__).parent / "schemas" / MANIFEST_SCHEMA_FILE, deposit / MANIFEST_SCHEMA_FILE
    )
    (deposit / "README.md").write_text(readme)
    tag = tag if commit else None
    prepared = prepare_upload(
        deposit,
        study_name=root.name,
        tag=tag,
        manifest=manifest,
        zenodo=json.loads((root / ZENODO).read_text()),
        warnings=warnings,
        kind=kind,
    )
    return FreezeResult(
        root,
        tag,
        commit,
        deposit,
        manifest,
        warnings,
        guide=prepared["guide"],
        upload=prepared["upload"],
    )


def freeze(
    root: str | Path,
    *,
    tag: str | None = None,
    publish: bool = True,
) -> FreezeResult:
    """Freeze the study in ``root`` for publication; see the module docstring.

    With ``publish=False``, as :func:`~polyzymd.analyses.project_freeze.freeze_project`
    freezes each study of a project, only the study's own files are written
    (``manifest.json``, ``md_checklist.yaml``, ``system_summary.csv``, and the
    engine inputs and final frames under ``deposit/``); the metadata,
    citation, commit, tag and upload are the project's.

    Raises
    ------
    ProtocolError
        Only when the study file or its metadata cannot be read, a condition
        config is outside the study folder and its project folder, or the tag
        already exists. Everything else is a warning in the result.
    """
    import yaml

    import polyzymd
    from polyzymd.analyses.study_file import load_study_file, outside_configs, portable
    from polyzymd.analyses.study_metadata import check_metadata

    protocol = load_study_file(root)
    root = protocol.root
    if publish and protocol.project is not None:
        raise ProtocolError(
            f"{root.name} is a study of the project {protocol.project.root}, whose "
            "project.yaml and shared analyses/ its results depend on.",
            hint=f"Freeze the whole project: polyzymd project freeze {protocol.project.root}",
        )
    for label, config in outside_configs(protocol).items():
        raise ProtocolError(
            f"The config of condition {label} is {config}, outside the study folder, "
            "so the deposit would not hold it.",
            hint=f"Remove {label} from conditions: in study.yaml, then copy it in with: "
            f'polyzymd study add-condition "{label}" --config {config}. Or move its folder '
            "under conditions/ and give the new path in study.yaml.",
        )
    meta, warnings = check_metadata(protocol.metadata)
    if publish:
        state, tag, found = _git_preflight(root, tag, "study")
        warnings += found
    else:
        # The project checks and publishes the metadata once, for every study,
        # and commits and tags them together; the study's files are those git
        # tracks in the project's repository.
        warnings, tag = [], None
        from polyzymd.analyses.study_git import git_state

        state = git_state(root)
    deposit = root / DEPOSIT
    for part in ("engine_inputs", "final_frames"):
        # A replicate no longer in the study must not stay in the deposit.
        shutil.rmtree(deposit / part, ignore_errors=True)
    deposit.mkdir(exist_ok=True)
    hashes = _Hashes()
    conditions, rows = _replicates(protocol, deposit, hashes, warnings)
    for run, why in stale_runs(protocol, conditions).items():
        warnings.append(f"run {run} may be stale: {'; '.join(why)}")
    warnings.extend(report_problems(protocol))
    warnings.extend(_outside_inputs(protocol))
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

    if publish:
        _write_gitignore(root)
    study_files = [p for p in _listed_files(root, state) if p not in GENERATED]
    deposit_root = protocol.project.root if protocol.project is not None else root
    configs = {path.resolve() for path in protocol.conditions.values()}
    # A project names the files it leaves out once, for all its studies.
    left_out = left_out_files(root, state) if publish else None
    if left_out:
        warnings.append(left_out)
    # The same warning for several conditions is one line naming them.
    warnings[:] = group_warnings(warnings, list(protocol.conditions))
    manifest: dict[str, Any] = {
        "$schema": MANIFEST_SCHEMA_FILE,
        "schema": MANIFEST_SCHEMA,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tag": tag,
        # The tag names the frozen commit, which holds this file and so
        # cannot be named in it; this is the commit before it.
        "git": {
            "parent_commit": state["commit"] if state else None,
            "inputs_uncommitted": state["inputs_uncommitted"] if state else None,
        },
        "study_file": {"path": STUDY_FILE, **hashes(protocol.path)},
        "versions": _versions(root),
        "equilibration": protocol.equilibration,
        "stride": protocol.stride,
        "metadata": meta,
        "conditions": conditions,
        "analyses": {
            run: {
                "analysis": entry.analysis,
                # A project's shared function lies outside the study: ../analyses/f.py.
                "function": f"{os.path.relpath(entry.function.file, root)}:{entry.function.qualname}"
                if entry.function
                else None,
                "settings": portable(
                    entry.settings or (entry.function.settings if entry.function else {}),
                    root,
                    protocol.project.root if protocol.project is not None else None,
                ),
                "selections": dict(entry.function.selections) if entry.function else {},
                "equilibration": protocol.window(run)[0],
                "until": protocol.window(run)[1],
                "stride": protocol.stride_of(run),
                "report": f"results/{run}/report.json",
            }
            for run, entry in protocol.analyses.items()
        },
        "trajectory_deposits": meta["related"]["trajectories"],
        # As deposited: a condition config without the machine paths it names.
        "files": {p: deposited_entry(root / p, deposit_root, configs) for p in study_files},
        "cite": {
            "polyzymd": __import__("polyzymd.citation", fromlist=["citation_line"]).citation_line()
        },
        "deposit_without_machine_paths": [
            f"{config}: machine paths removed; its entry in files is the deposited copy"
            for config in _condition_configs(protocol)
            if (root / config).is_file()
            and deposited_entry(root / config, deposit_root, configs) != hashes(root / config)
        ],
        "warnings": warnings,
    }
    (root / CHECKLIST).write_text(
        "# The Communications Biology MD checklist, filled by polyzymd study freeze. Informational.\n"
        + yaml.safe_dump(_checklist(protocol, manifest, meta), sort_keys=False)
    )
    if publish:
        _write_citation(
            root,
            meta,
            version=version,
            released=released,
            commit=None,
            method=_method(protocol),
        )
    # The files freeze writes are deposited too; the manifest cannot list itself.
    for name in (CHECKLIST, SUMMARY, *((CITATION, ZENODO) if publish else ())):
        if (root / name).is_file():
            manifest["files"][name] = hashes(root / name)
    (root / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    if not publish:
        return FreezeResult(root, None, None, deposit, manifest, warnings)
    commit = None
    if state and tag:
        paths = [p for p in (*GENERATED, ".gitignore", "results") if (root / p).exists()]
        commit = _commit_and_tag(root, paths, tag, "study", warnings)
        if commit is None:
            _drop_tag(root, manifest, meta, released, _method(protocol))
    generated = [p for p in GENERATED if (root / p).exists()]
    files = sorted({*study_files, *generated})
    _copy_frozen_folder(root, deposit, tag, commit, files, _condition_configs(protocol))
    from polyzymd.analyses.study_upload_guide import deposit_readme

    # The deposit's README describes the study from its metadata; the study's
    # own README.md stays as written, inside the study.
    readme = deposit_readme(
        study_name=root.name,
        tag=tag if commit else None,
        meta=meta,
        analyses=protocol.analyses,
        root=root,
    )
    result = _finish_deposit(root, deposit, tag, commit, manifest, readme, warnings)
    result.git_failed = bool(state and tag and commit is None)
    return result
