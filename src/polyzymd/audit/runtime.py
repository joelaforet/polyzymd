"""At run time: is the System about to be integrated the one the build audited?

The build audit vouches for ``system.xml``. A run adds what a run should (a
barostat, position restraints for a stage) and nothing else; anything more is
a different system from the one audited, as when a production run once
recomputed charges. Every force-field term is compared; run controls and
restraint forces are set aside, as the audit sets them aside.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from polyzymd.audit.terms import TERM_KINDS, SystemIndex, compare

#: Name of the per-segment file the comparison is written to.
SYSTEM_CHECK_SUFFIX = "_system_check.json"


def _restraints(system: Any) -> dict[str, Any]:
    """Every restraint-type force, by name: its energy and each bond's atoms and parameters."""
    from polyzymd.audit.terms import RESTRAINT_FORCES

    found = {}
    for number, force in enumerate(system.getForces()):
        if force.__class__.__name__ not in RESTRAINT_FORCES:
            continue
        name = force.getName() or f"{force.__class__.__name__} {number}"
        bonds = []
        if hasattr(force, "getNumBonds"):
            for k in range(force.getNumBonds()):
                *atoms, values = force.getBondParameters(k)
                bonds.append((tuple(atoms), tuple(round(float(v), 9) for v in values)))
        energy = force.getEnergyFunction() if hasattr(force, "getEnergyFunction") else ""
        found[name] = (energy, sorted(bonds))
    return found


def compare_with_build(system: Any, build_system_xml: Path) -> dict[str, Any]:
    """Compare ``system`` with the audited ``system.xml`` of the build.

    The build record next to ``system.xml`` must vouch for it: a passing audit
    and the file's recorded hash. Then every force-field term, the restraints
    and the nonbonded settings are compared; run controls (barostat, motion
    remover) are not. Returns a plain dict: ``compared`` (whether the
    comparison ran), ``matches_build``, ``detail`` and the differences.
    """
    import json

    import openmm

    from polyzymd.analyses.shared.file_hashes import file_sha256
    from polyzymd.audit.realized import _nonbonded_settings
    from polyzymd.simulation.artifact_integrity import build_record_path

    build_system_xml = Path(build_system_xml)
    record: dict[str, Any] = {"build_system": str(build_system_xml), "compared": False}
    if not build_system_xml.is_file():
        record.update(matches_build=False, detail=f"{build_system_xml} is missing")
        return record
    sha = file_sha256(build_system_xml, use_cache=False)
    record["build_sha256"] = sha
    path = build_record_path(build_system_xml.parent)
    try:
        build = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        build = {}
    audit = build.get("audit") or {}
    recorded = (build.get("artifacts") or {}).get("system.xml", {}).get("sha256")
    if not audit.get("passed") or recorded != sha:
        reason = (
            "its build record has no passing physics audit"
            if not audit.get("passed")
            else "it changed since the build was audited"
        )
        record.update(
            matches_build=False,
            detail=f"{build_system_xml} cannot vouch for this run: {reason} ({path.name})",
        )
        return record
    built = openmm.XmlSerializer.deserialize(build_system_xml.read_text())
    if built.getNumParticles() != system.getNumParticles():
        record.update(
            compared=True,
            matches_build=False,
            detail=f"{system.getNumParticles()} particles, the build has {built.getNumParticles()}",
        )
        return record
    atoms = range(system.getNumParticles())
    running = SystemIndex.from_system(system).extract(atoms)
    audited = SystemIndex.from_system(built).extract(atoms)
    differing: dict[str, Any] = {}
    for kind in TERM_KINDS:
        result = compare(kind, running[kind], audited[kind])
        if not result.ok:
            differing[kind] = result.total - result.matched
    now, then = _restraints(system), _restraints(built)
    changed = sorted(n for n in set(now) | set(then) if now.get(n) != then.get(n))
    if changed:
        differing["restraints"] = changed
    settings_now, settings_then = _nonbonded_settings(system), _nonbonded_settings(built)
    settings = sorted(k for k in settings_then if settings_now.get(k) != settings_then.get(k))
    if settings:
        differing["nonbonded_settings"] = settings
    record["compared"] = True
    record["differing_terms"] = differing
    record["matches_build"] = not differing
    record["detail"] = (
        "every force-field term, restraint and nonbonded setting matches the audited build"
        if not differing
        else "differs from the audited build: "
        + ", ".join(
            f"{v} {k}" if isinstance(v, int) else f"{k} {', '.join(map(str, v))}"
            for k, v in differing.items()
        )
    )
    return record


class SystemDivergedError(RuntimeError):
    """Raised when a segment would integrate a System other than the audited build."""


def record_system_check(
    system: Any, working_dir: Path, phase_dir: Path, phase_name: str, segment_index: int
) -> dict[str, Any]:
    """Compare the segment's System with the audited build and record the result.

    Writes ``<phase>_system_check.json`` beside the segment and the same record
    on the segment in ``progress.json`` (shown by ``polyzymd status``). A
    System that differs from the audited build stops the segment before its
    first step (:class:`SystemDivergedError`). When no comparison could run
    (no passing build record vouches for ``system.xml``, as for a run started
    before 1.3.0rc5), the record says so and the run goes on.
    """
    import json
    import logging

    from polyzymd.simulation.progress import update_segment_fields

    logger = logging.getLogger(__name__)
    try:
        record = compare_with_build(system, Path(working_dir) / "system.xml")
    except Exception as error:  # a comparison that cannot run is not a match
        record = {
            "matches_build": False,
            "detail": f"could not compare with the build: {type(error).__name__}: {error}",
        }
    (Path(phase_dir) / f"{phase_name}{SYSTEM_CHECK_SUFFIX}").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n"
    )
    if record["matches_build"]:
        logger.info(f"{phase_name}: the system matches the audited build")
    else:
        logger.error(
            f"{phase_name}: THE SYSTEM BEING RUN IS NOT THE AUDITED BUILD: {record['detail']}"
        )
    update_segment_fields(working_dir, segment_index, system_check=record)
    if record.get("compared") and not record["matches_build"]:
        raise SystemDivergedError(
            f"SYSTEM DIVERGED FROM THE AUDITED BUILD: {phase_name} would integrate a System "
            f"that {record['detail']}; stopped before its first step (the record is in "
            f"{Path(phase_dir) / (phase_name + SYSTEM_CHECK_SUFFIX)}).\n"
            "fix: find what changed the System since the build (the config, the OpenMM or "
            "PolyzyMD version, a hand-edited file) and build the replicate again, or run it "
            "from the audited build unchanged. Do not resubmit it as it is."
        )
    return record
