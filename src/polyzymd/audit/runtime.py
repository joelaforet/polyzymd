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


def compare_with_build(system: Any, build_system_xml: Path) -> dict[str, Any]:
    """Compare ``system`` with the audited ``system.xml`` of the build.

    Returns a plain dict: ``matches_build`` (bool), ``detail``, the build file's
    SHA-256, and per kind the number of differing terms.
    """
    import openmm

    from polyzymd.analyses.shared.file_hashes import file_sha256

    build_system_xml = Path(build_system_xml)
    record: dict[str, Any] = {"build_system": str(build_system_xml)}
    if not build_system_xml.is_file():
        record.update(matches_build=False, detail=f"{build_system_xml} is missing")
        return record
    record["build_sha256"] = file_sha256(build_system_xml, use_cache=False)
    built = openmm.XmlSerializer.deserialize(build_system_xml.read_text())
    if built.getNumParticles() != system.getNumParticles():
        record.update(
            matches_build=False,
            detail=f"{system.getNumParticles()} particles, the build has {built.getNumParticles()}",
        )
        return record
    atoms = range(system.getNumParticles())
    running = SystemIndex.from_system(system).extract(atoms)
    audited = SystemIndex.from_system(built).extract(atoms)
    differing = {}
    for kind in TERM_KINDS:
        result = compare(kind, running[kind], audited[kind])
        if not result.ok:
            differing[kind] = result.total - result.matched
    record["differing_terms"] = differing
    record["matches_build"] = not differing
    record["detail"] = (
        "every force-field term matches the audited build"
        if not differing
        else "differs from the audited build: "
        + ", ".join(f"{n} {kind}" for kind, n in differing.items())
    )
    return record


def record_system_check(
    system: Any, working_dir: Path, phase_dir: Path, phase_name: str, segment_index: int
) -> dict[str, Any]:
    """Compare the segment's System with the audited build and record the result.

    Writes ``<phase>_system_check.json`` beside the segment and the same record
    on the segment in ``progress.json`` (shown by ``polyzymd status``). A
    divergence is logged as an error; it does not stop the run.
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
    return record
