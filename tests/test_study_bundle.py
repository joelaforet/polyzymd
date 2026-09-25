"""Packaging a study for publication and checking the package someone received."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.analyses.identity import file_fingerprint
from polyzymd.study_bundle import MANIFEST_NAME, export_study, plan_export, verify_study


@pytest.fixture
def study(tmp_path: Path) -> Path:
    """A study with one referenced condition, one rerun nobody lists, and a trajectory."""
    root = tmp_path / "thermal"
    (root / "analyses").mkdir(parents=True)
    (root / "study.yaml").write_text("name: thermal\n")
    (root / "analyses" / "lid_opening.py").write_text("# analysis\n")
    for condition in ("noPoly_343K", "noPoly_343K_REDO"):
        directory = root / "conditions" / condition
        (directory / "structures").mkdir(parents=True)
        (directory / "config.yaml").write_text(f"name: {condition}\n")
        (directory / "structures" / "enzyme.pdb").write_text("ATOM\n")
    run_dir = root / "conditions" / "noPoly_343K" / "run_1" / "production_0"
    run_dir.mkdir(parents=True)
    trajectory = run_dir / "prod.dcd"
    trajectory.write_bytes(b"DCD" * 100)
    (root / "conditions" / "noPoly_343K" / "slurm_logs").mkdir()
    (root / "conditions" / "noPoly_343K" / "slurm_logs" / "job.out").write_text("log\n")

    comparison = root / "comparisons" / "calb_343K"
    comparison.mkdir(parents=True)
    (comparison / "comparison.yaml").write_text(
        "name: CALB 343 K\n"
        "conditions:\n"
        "  - label: No Polymer\n"
        "    config: ../../conditions/noPoly_343K/config.yaml\n"
        "    replicates: [1]\n"
    )
    result = comparison / "analysis" / "No_Polymer" / "rg" / "run_1" / "result.json"
    result.parent.mkdir(parents=True)
    result.write_text(
        json.dumps(
            {
                "condition_label": "No Polymer",
                "replicate": 1,
                "provenance": {
                    "identity": {
                        "inputs": [
                            {
                                "path": str(trajectory),
                                "relative_path": "production_0/prod.dcd",
                                "size_bytes": trajectory.stat().st_size,
                                "fingerprint": file_fingerprint(trajectory),
                            }
                        ]
                    }
                },
            }
        )
    )
    return root


def _working_dir(root: Path):
    return lambda config, replicate: config.parent / f"run_{replicate}"


def test_the_package_holds_what_the_comparisons_reference(study: Path) -> None:
    plan = plan_export(study)
    packaged = {path.relative_to(study).as_posix() for path in plan.files}

    assert "study.yaml" in packaged
    assert "analyses/lid_opening.py" in packaged
    assert "conditions/noPoly_343K/config.yaml" in packaged
    assert "conditions/noPoly_343K/structures/enzyme.pdb" in packaged
    assert "comparisons/calb_343K/comparison.yaml" in packaged
    assert not any("REDO" in path for path in packaged)
    assert not any(path.endswith(".dcd") for path in packaged)
    assert not any("slurm_logs" in path for path in packaged)
    assert [path.name for path in plan.unreferenced] == ["noPoly_343K_REDO"]


def test_the_manifest_lists_the_trajectories_to_archive(study: Path, tmp_path: Path) -> None:
    export_study(study, tmp_path / "thermal.zip")

    with zipfile.ZipFile(tmp_path / "thermal.zip") as archive:
        manifest = json.loads(archive.read(f"thermal/{MANIFEST_NAME}"))

    assert manifest["study"] == "thermal"
    assert manifest["trajectories"] == [
        {
            "condition": "No Polymer",
            "config": "conditions/noPoly_343K/config.yaml",
            "replicate": 1,
            "relative_path": "production_0/prod.dcd",
            "size_bytes": 300,
            "fingerprint": file_fingerprint(
                study / "conditions" / "noPoly_343K" / "run_1" / "production_0" / "prod.dcd"
            ),
        }
    ]
    assert manifest["left_out_conditions"] == ["conditions/noPoly_343K_REDO"]


def test_an_unpacked_package_verifies(study: Path, tmp_path: Path) -> None:
    export_study(study, tmp_path / "thermal.zip")
    with zipfile.ZipFile(tmp_path / "thermal.zip") as archive:
        archive.extractall(tmp_path / "received")
    received = tmp_path / "received" / "thermal"

    report = verify_study(received, working_dir=_working_dir(received))

    assert report.ok
    assert report.checked_files > 5
    assert len(report.trajectories_absent) == 1


def test_an_edited_or_missing_file_fails_verification(study: Path, tmp_path: Path) -> None:
    export_study(study, tmp_path / "thermal.zip")
    with zipfile.ZipFile(tmp_path / "thermal.zip") as archive:
        archive.extractall(tmp_path / "received")
    received = tmp_path / "received" / "thermal"
    (received / "analyses" / "lid_opening.py").write_text("# edited\n")
    (received / "study.yaml").unlink()

    report = verify_study(received, working_dir=_working_dir(received))

    assert not report.ok
    assert report.changed_files == ["analyses/lid_opening.py"]
    assert report.missing_files == ["study.yaml"]


def test_downloaded_trajectories_are_checked_by_content(study: Path, tmp_path: Path) -> None:
    export_study(study, tmp_path / "thermal.zip")
    with zipfile.ZipFile(tmp_path / "thermal.zip") as archive:
        archive.extractall(tmp_path / "received")
    received = tmp_path / "received" / "thermal"
    downloaded = received / "conditions" / "noPoly_343K" / "run_1" / "production_0"
    downloaded.mkdir(parents=True)

    (downloaded / "prod.dcd").write_bytes(b"DCD" * 100)
    good = verify_study(received, working_dir=_working_dir(received))
    (downloaded / "prod.dcd").write_bytes(b"XYZ" * 100)
    bad = verify_study(received, working_dir=_working_dir(received))

    assert good.ok and len(good.trajectories_ok) == 1
    assert not bad.ok and len(bad.trajectories_changed) == 1


def test_export_leaves_out_its_own_zip(study: Path) -> None:
    export_study(study, study / "thermal.zip")
    export_study(study, study / "thermal.zip")

    with zipfile.ZipFile(study / "thermal.zip") as archive:
        assert "thermal/thermal.zip" not in archive.namelist()


def test_the_export_and_verify_commands(
    study: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from polyzymd.cli.main import cli

    monkeypatch.chdir(study / "comparisons")
    exported = CliRunner().invoke(cli, ["study", "export", "-o", str(tmp_path / "t.zip")])
    with zipfile.ZipFile(tmp_path / "t.zip") as archive:
        archive.extractall(tmp_path / "received")
    verified = CliRunner().invoke(cli, ["study", "verify", str(tmp_path / "received" / "thermal")])

    assert exported.exit_code == 0, exported.output
    assert "1 trajectory file(s) to archive separately" in exported.output
    assert "Left out conditions/noPoly_343K_REDO" in exported.output
    assert verified.exit_code == 0, verified.output
    assert "not downloaded" in verified.output
