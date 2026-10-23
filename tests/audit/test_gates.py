"""The gates: a failed audit stops the build, and runs refuse files no passing audit vouches for."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

openmm = pytest.importorskip("openmm")

from polyzymd.audit.findings import Finding  # noqa: E402
from polyzymd.audit.report import AuditReport  # noqa: E402
from polyzymd.simulation import artifact_integrity  # noqa: E402
from polyzymd.simulation.artifact_integrity import (  # noqa: E402
    MANIFEST_NAME,
    AuditFailedError,
    AuditGateError,
    config_hash,
    require_audited_build,
)


class _Config:
    """Just enough config for config_hash."""

    def __init__(self, tag: str = "a") -> None:
        self.tag = tag

    def model_dump(self, **kwargs):
        return {"tag": self.tag}


def _record(folder: Path, config, *, passed=True, engine="openmm", schema=2, audit=True):
    (folder / "system.xml").write_text("<System/>")
    from polyzymd.analyses.shared.file_hashes import file_sha256

    record = {
        "schema_version": schema,
        "engine": engine,
        "config_hash": config_hash(config),
        "artifacts": {
            "system.xml": {
                "path": "system.xml",
                "sha256": file_sha256(folder / "system.xml", use_cache=False),
            }
        },
        "audit": (
            {
                "passed": passed,
                "findings": (
                    [
                        {
                            "check": "parameters.bonded",
                            "component": "protein",
                            "status": "FAIL",
                            "detail": "Sage on protein",
                        }
                    ]
                    if not passed
                    else []
                ),
            }
            if audit
            else None
        ),
    }
    (folder / MANIFEST_NAME).write_text(json.dumps(record))
    return record


def test_a_passing_audited_build_is_accepted(tmp_path) -> None:
    config = _Config()
    _record(tmp_path, config)
    assert require_audited_build(tmp_path, config, "openmm")["audit"]["passed"]


@pytest.mark.parametrize(
    "case, match",
    [
        ("missing", "has no build record"),
        ("legacy", "records no physics audit"),
        ("no_audit", "records no physics audit"),
        ("engine", "records a gromacs build"),
        ("config", "config changed"),
        ("file", "missing or changed"),
        ("failed", "physics audit failed"),
    ],
)
def test_a_run_refuses_files_no_passing_audit_vouches_for(tmp_path, case, match) -> None:
    config = _Config()
    if case != "missing":
        _record(
            tmp_path,
            config,
            schema=1 if case == "legacy" else 2,
            audit=case != "no_audit",
            engine="gromacs" if case == "engine" else "openmm",
            passed=case != "failed",
        )
    if case == "config":
        config = _Config("b")
    if case == "file":
        (tmp_path / "system.xml").write_text("<System changed/>")
    with pytest.raises(AuditGateError, match=match) as raised:
        require_audited_build(tmp_path, config, "openmm")
    assert "fix:" in str(raised.value)
    if case == "failed":
        assert "parameters.bonded (protein): Sage on protein" in str(raised.value)
        assert "modeling_exceptions" in str(raised.value)


def _tiny_bundle_args(tmp_path):
    from openmm import app, unit

    topology = app.Topology()
    residue = topology.addResidue("X", topology.addChain())
    topology.addAtom("C", app.element.carbon, residue)
    system = openmm.System()
    system.addParticle(12.0)
    force = openmm.NonbondedForce()
    force.addParticle(0.0, 0.3, 0.1)
    system.addForce(force)
    return topology, system, [openmm.Vec3(0, 0, 0)] * unit.nanometer


class _Phases:
    class production:
        time_step = 2.0


class _BuildConfig(_Config):
    simulation_phases = _Phases


class _Component:
    def to_dict(self):
        return {
            "name": "protein",
            "type": "protein",
            "segid": "PROT",
            "molecules": [],
            "domains": [],
        }


def test_a_failed_audit_commits_nothing_and_keeps_the_record(tmp_path, monkeypatch) -> None:
    failing = AuditReport([Finding("parameters.bonded", "protein", "FAIL", "Sage on protein")])
    monkeypatch.setattr(artifact_integrity, "_audit", lambda *args: failing)
    topology, system, positions = _tiny_bundle_args(tmp_path)
    with pytest.raises(AuditFailedError, match="Sage on protein"):
        artifact_integrity.publish_build_bundle(
            tmp_path, topology, system, positions, _BuildConfig(), components=[_Component()]
        )
    assert not (tmp_path / "system.xml").exists()
    assert not (tmp_path / MANIFEST_NAME).exists()
    failed = json.loads((tmp_path / artifact_integrity.FAILED_RECORD).read_text())
    assert failed["audit"]["passed"] is False
    assert (tmp_path / "build_provenance.failed.md").is_file()


def _system(charge: float):
    system = openmm.System()
    force = openmm.NonbondedForce()
    for q in (charge, -charge):
        system.addParticle(12.0)
        force.addParticle(q, 0.3, 0.1)
    system.addForce(force)
    return system


def test_a_run_whose_system_differs_from_the_build_is_recorded_and_shown(tmp_path) -> None:
    """The LipA case: production ran with charges other than the audited build's."""
    from polyzymd.audit.runtime import compare_with_build, record_system_check
    from polyzymd.cli.status_report import (
        ReplicateReport,
        _diverged_segments,
        render_replicate_line,
    )

    build = tmp_path / "system.xml"
    build.write_text(openmm.XmlSerializer.serialize(_system(0.5)))
    same = _system(0.5)
    same.addForce(openmm.MonteCarloBarostat(1.0, 300.0))  # a run control, as production adds
    assert compare_with_build(same, build)["matches_build"]
    record = compare_with_build(_system(0.4), build)
    assert not record["matches_build"] and record["differing_terms"] == {"charge": 2}

    phase = tmp_path / "production_0"
    phase.mkdir()
    written = record_system_check(_system(0.4), tmp_path, phase, "production_0", 0)
    assert json.loads((phase / "production_0_system_check.json").read_text()) == written

    class _Segment:
        index = 0
        system_check = written

    class _Progress:
        segments = [_Segment()]

    diverged = _diverged_segments(_Progress())
    line = render_replicate_line(
        ReplicateReport(1, str(tmp_path), "running", 1.0, 10.0, 0.1, "running", diverged=diverged),
        4,
    )
    assert "SYSTEM DIVERGED FROM THE AUDITED BUILD" in line and "charge" in line


def test_progress_carries_the_segments_system_check(tmp_path) -> None:
    """The check is written when a segment starts; progress.json picks it up on save."""
    from polyzymd.simulation.progress import SegmentRecord, SimulationProgress, save_progress

    phase = tmp_path / "production_0"
    phase.mkdir()
    (phase / "production_0_system_check.json").write_text(
        json.dumps({"matches_build": False, "detail": "differs: 2 charge"})
    )
    progress = SimulationProgress(
        total_steps_requested=10, total_samples_requested=1, segments=[SegmentRecord(index=0)]
    )
    save_progress(tmp_path, progress)
    saved = json.loads((tmp_path / "progress.json").read_text())
    assert saved["segments"][0]["system_check"]["matches_build"] is False


def test_gromacs_submission_needs_an_audited_build_until_the_run_starts(tmp_path) -> None:
    from polyzymd.engines.base import EngineSubmitRequest
    from polyzymd.engines.gromacs.engine import GromacsEngine

    config = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "config_gromacs.yaml"
    gromacs_dir = tmp_path / "run_1" / "gromacs"
    gromacs_dir.mkdir(parents=True)
    request = EngineSubmitRequest(replicate=1, config_path=config, working_dir=gromacs_dir)
    with pytest.raises(AuditGateError, match="has no build record"):
        GromacsEngine._require_audited_build(request)
    (gromacs_dir / "em.gro").write_text("")  # started: continues from its own files
    GromacsEngine._require_audited_build(request)
