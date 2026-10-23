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


def _vouch(folder: Path, passed: bool = True) -> None:
    """A build record whose audit vouches (or not) for folder/system.xml."""
    from polyzymd.analyses.shared.file_hashes import file_sha256

    sha = file_sha256(folder / "system.xml", use_cache=False)
    (folder / MANIFEST_NAME).write_text(
        json.dumps(
            {
                "schema_version": 2,
                "artifacts": {"system.xml": {"path": "system.xml", "sha256": sha}},
                "audit": {"passed": passed, "findings": []},
            }
        )
    )


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
    _vouch(tmp_path)
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


def test_the_runtime_check_compares_restraints_and_needs_a_vouching_record(tmp_path) -> None:
    from polyzymd.audit.runtime import compare_with_build

    def restrained(k):
        system = _system(0.5)
        force = openmm.CustomBondForce("0.5 * k * (r - r0)^2")
        force.addPerBondParameter("r0")
        force.addPerBondParameter("k")
        force.addBond(0, 1, [0.3, k])
        force.setName("polyzymd_restraint:tether")
        system.addForce(force)
        return system

    build = tmp_path / "system.xml"
    build.write_text(openmm.XmlSerializer.serialize(restrained(1000.0)))
    _vouch(tmp_path)
    assert compare_with_build(restrained(1000.0), build)["matches_build"]
    changed = compare_with_build(restrained(2000.0), build)
    assert not changed["matches_build"] and "polyzymd_restraint:tether" in changed["detail"]
    _vouch(tmp_path, passed=False)
    unvouched = compare_with_build(restrained(1000.0), build)
    assert (unvouched["compared"], unvouched["matches_build"]) == (False, False)
    assert "no passing physics audit" in unvouched["detail"]


def test_a_started_replicate_is_not_held_to_the_audit_it_predates(tmp_path) -> None:
    from polyzymd.simulation.artifact_integrity import run_started

    assert not run_started(tmp_path)
    (tmp_path / "production_0").mkdir()
    assert run_started(tmp_path)


def test_the_gromacs_job_refuses_files_no_audit_vouched_for() -> None:
    """The job script checks audited_files.sha256 before any grompp."""
    template = (
        Path(__file__).resolve().parents[2]
        / "src/polyzymd/engines/gromacs/templates/slurm_self_resubmitting.sh.jinja"
    ).read_text()
    check = template.index("sha256sum --check --quiet audited_files.sha256")
    assert check < template.index("$GMX grompp")
    assert "no passing physics audit vouches for these files" in template


def _hash_check_block(template_name: str) -> str:
    text = (
        Path(__file__).resolve().parents[2]
        / "src/polyzymd/engines/gromacs/templates"
        / template_name
    ).read_text()
    start = text.index("# Run only files a passing physics audit vouched for")
    end = text.index("\nfi\n", start) + len("\nfi\n")
    return text[start:end]


def test_both_gromacs_scripts_check_the_audited_files_the_same_way() -> None:
    assert _hash_check_block("run_gromacs.sh.jinja") == _hash_check_block(
        "slurm_self_resubmitting.sh.jinja"
    )


@pytest.mark.parametrize(
    "case, runs",
    [
        ("not started, no hashes", False),
        ("started before rc.5", True),
        ("edited after the audit", False),
        ("as audited", True),
    ],
)
def test_the_gromacs_hash_check_refuses_unvouched_files_but_not_old_runs(
    tmp_path, case, runs
) -> None:
    import subprocess

    (tmp_path / "prod.mdp").write_text("dt = 0.002\n")
    if case == "started before rc.5":
        (tmp_path / "em.gro").write_text("")
    if case in ("edited after the audit", "as audited"):
        subprocess.run(
            "sha256sum prod.mdp > audited_files.sha256", shell=True, cwd=tmp_path, check=True
        )
    if case == "edited after the audit":
        (tmp_path / "prod.mdp").write_text("dt = 0.004\n")
    script = _hash_check_block("run_gromacs.sh.jinja") + "echo RAN\n"
    result = subprocess.run(["bash", "-c", script], cwd=tmp_path, capture_output=True, text=True)
    assert ("RAN" in result.stdout) is runs, result.stdout + result.stderr


def test_a_config_without_nagl_names_no_nagl_model() -> None:
    from polyzymd.config.loader import load_config

    root = Path(__file__).resolve().parents[2]
    assert load_config(root / "examples/quickstart/config.yaml").nagl_model is None
    with_nagl = load_config(root / "tests/config/data/rc5_configs/enzyme_cosolvent.yaml")
    models = {
        c.charges.model for c in with_nagl.solvent.co_solvents if c.charges.method == "nagl"
    } | {s.charges.model for s in with_nagl.solutes if s.charges.method == "nagl"}
    assert with_nagl.nagl_model == (models.pop() if models else None)
