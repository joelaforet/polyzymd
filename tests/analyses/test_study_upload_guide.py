"""deposit/upload/, deposit/trajectories.csv and deposit/UPLOAD.md, written by study freeze."""

from __future__ import annotations

import csv
import zipfile
from pathlib import Path

import pytest

from polyzymd.analyses.study_upload_guide import (
    FILE_BYTES,
    RECORD_BYTES,
    RECORD_FILES,
    prepare_upload,
    trajectory_batches,
)

GB = 1_000_000_000


def _manifest(conditions: dict[str, dict[str, list[int]]]) -> dict:
    """A manifest whose replicates hold files of the given sizes."""
    return {
        "metadata": {"doi": None},
        "conditions": {
            label: {
                "replicates": {
                    str(index): {
                        "files": [
                            {"path": f"run{index}/part{k}.dcd", "size": size, "sha256": "0" * 64}
                            for k, size in enumerate(sizes)
                        ]
                    }
                    for index, sizes in replicates.items()
                }
            }
            for label, replicates in conditions.items()
        },
    }


ZENODO = {
    "title": "T",
    "publication_date": "2026-10-01",
    "creators": [{"name": "Lovelace, Ada", "orcid": "0000-0002-1825-0097"}],
    "description": "D",
    "license": "cc-by-4.0",
    "keywords": ["md"],
    "version": "study-v1",
    "method": "M",
    "related_identifiers": [
        {
            "identifier": "https://github.com/joelaforet/polyzymd",
            "relation": "requires",
            "resource_type": "software",
        }
    ],
}


class TestBatches:
    def test_replicates_share_a_record_until_it_is_full(self) -> None:
        batches, notes = trajectory_batches(
            _manifest({"A": {1: [20 * GB], 2: [20 * GB], 3: [20 * GB]}})
        )
        assert [b.replicates for b in batches] == [("1", "2"), ("3",)]
        assert all(b.size <= RECORD_BYTES for b in batches) and notes == []

    def test_conditions_never_share_a_batch(self) -> None:
        batches, _ = trajectory_batches(_manifest({"A": {1: [GB]}, "B": {1: [GB]}}))
        assert [b.condition for b in batches] == ["A", "B"]

    def test_file_count_limit(self) -> None:
        batches, _ = trajectory_batches(_manifest({"A": {1: [1] * 60, 2: [1] * 60}}))
        assert [len(b.files) for b in batches] == [60, 60]
        assert all(len(b.files) <= RECORD_FILES for b in batches)

    def test_oversized_file_and_replicate_are_named(self) -> None:
        _, notes = trajectory_batches(_manifest({"A": {1: [FILE_BYTES + 1, 30 * GB]}}))
        assert any("larger than Zenodo's" in n for n in notes)
        assert any("quota" in n for n in notes)


class TestPrepare:
    def _deposit(self, tmp_path: Path) -> Path:
        deposit = tmp_path / "deposit"
        for name in ("README.md", "CITATION.cff", "manifest.json"):
            (deposit / name).parent.mkdir(parents=True, exist_ok=True)
            (deposit / name).write_text(name)
        (deposit / "study" / "results").mkdir(parents=True)
        (deposit / "study" / "study.yaml").write_text("x")
        (deposit / "engine_inputs" / "a").mkdir(parents=True)
        (deposit / "engine_inputs" / "a" / "system.xml.gz").write_bytes(b"x")
        return deposit

    def test_upload_folder(self, tmp_path: Path) -> None:
        deposit = self._deposit(tmp_path)
        out = prepare_upload(
            deposit,
            study_name="s",
            tag="study-v1",
            manifest=_manifest({"A": {1: [GB]}}),
            zenodo=ZENODO,
            warnings=[],
        )
        names = sorted(p.name for p in out["upload"].iterdir())
        assert names == [
            "CITATION.cff",
            "README.md",
            "engine_inputs.zip",
            "manifest.json",
            "s-study-v1.zip",
        ]
        assert len(names) <= RECORD_FILES
        with zipfile.ZipFile(out["upload"] / "s-study-v1.zip") as archive:
            assert "study/study.yaml" in archive.namelist()

    def test_trajectories_csv(self, tmp_path: Path) -> None:
        out = prepare_upload(
            self._deposit(tmp_path),
            study_name="s",
            tag="v",
            manifest=_manifest({"A": {1: [GB, 2]}}),
            zenodo=ZENODO,
            warnings=[],
        )
        rows = list(csv.DictReader(out["trajectories"].open()))
        assert [r["path"] for r in rows] == ["run1/part0.dcd", "run1/part1.dcd"]
        assert (
            rows[0]["batch"] == "1"
            and rows[0]["size_bytes"] == str(GB)
            and len(rows[0]["sha256"]) == 64
        )

    def test_guide_without_a_doi(self, tmp_path: Path) -> None:
        out = prepare_upload(
            self._deposit(tmp_path),
            study_name="s",
            tag="study-v1",
            manifest=_manifest({"A": {1: [GB]}}),
            zenodo=ZENODO,
            warnings=["metadata.purpose is missing"],
        )
        text = out["guide"].read_text()
        assert "Get a DOI now!" in text and "metadata.purpose is missing" in text
        assert "| Creators | Lovelace, Ada (ORCID 0000-0002-1825-0097) |" in text
        assert "requires https://github.com/joelaforet/polyzymd" in text
        assert "| 1 | A | 1 | 1 | 1.00 GB |" in text
        assert (
            "PolyzyMD does not upload trajectories" in text
            and "uploading and publishing are yours" in text
        )

    def test_guide_with_a_doi(self, tmp_path: Path) -> None:
        manifest = _manifest({"A": {1: [GB]}})
        manifest["metadata"]["doi"] = "10.5281/zenodo.1234567"
        out = prepare_upload(
            self._deposit(tmp_path),
            study_name="s",
            tag="v",
            manifest=manifest,
            zenodo=ZENODO,
            warnings=[],
        )
        text = out["guide"].read_text()
        assert "already set: `10.5281/zenodo.1234567`" in text and "- none" in text

    def test_no_trajectories(self, tmp_path: Path) -> None:
        out = prepare_upload(
            self._deposit(tmp_path),
            study_name="s",
            tag="v",
            manifest={"metadata": {}, "conditions": {}},
            zenodo=ZENODO,
            warnings=[],
        )
        assert "no trajectories were on this machine" in out["guide"].read_text()


@pytest.mark.parametrize(
    ("size", "text"),
    [(999, "999 B"), (1500, "1.50 kB"), (2_500_000, "2.50 MB"), (3 * GB, "3.00 GB")],
)
def test_sizes(size: int, text: str) -> None:
    from polyzymd.analyses.study_upload_guide import _gb

    assert _gb(size) == text


@pytest.mark.parametrize("project", [True, False])
def test_deposit_readme_says_where_to_locate_and_how_to_analyze_one_study(
    tmp_path: Path, project: bool
) -> None:
    """The README says which folder runs `study locate`, and a project's
    README gives the per-study analyze for a reader who downloaded one study."""
    from polyzymd.analyses.study_upload_guide import deposit_readme

    text = deposit_readme(
        study_name="paper", tag="v1", meta={}, analyses={}, root=tmp_path, project=project
    )
    assert "the unzipped `study/` folder" in text
    assert "DOWNLOAD_DIR is the folder that holds" in text
    if project:
        assert "`polyzymd study locate DOWNLOAD_DIR --verify --study <study>`" in text
        assert "`polyzymd analyze --project . --recompute`" in text
        assert "`polyzymd analyze --study <study> --recompute`" in text
    else:
        assert "`polyzymd study locate DOWNLOAD_DIR --verify`" in text
        assert "`polyzymd analyze --study study.yaml --recompute`" in text
    assert "Without `--recompute`, analyze reads the stored results" in text
