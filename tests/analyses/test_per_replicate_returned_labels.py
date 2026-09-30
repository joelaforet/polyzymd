"""Tests for ``Study.per_replicate(..., labels="returned")``.

Each replicate is an OpenMM run directory of the four-atom cross of
tests/_support/analysis_testkit.py, scaled by a number that differs between
replicates and conditions. :func:`_measure` reads that number from the
first production frame and returns the labels and values that :data:`TABLE`
gives it, so the labels differ between replicates. :data:`CALLS` counts the
calls, to tell a computed replicate from a reused one.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

#: Scale of the cross -> (labels, values) that replicate returns.
TABLE = {
    1.0: (["12-SBM", "45-EGM"], [0.5, 0.25]),
    2.0: (["45-EGM", "7-SBM"], [0.75, 0.1]),
    3.0: ([], []),
    4.0: (["12-SBM"], [0.9]),
    5.0: (["7-SBM", "12-SBM"], [0.2, 0.3]),
    6.0: ([np.int64(3), np.int64(1)], [0.4, 0.6]),
}
#: Condition -> scale of each replicate.
SCALES = {"A": (1.0, 2.0, 3.0), "B": (4.0, 5.0, 1.0)}
CALLS: list[float] = []


def _measure(atoms, frames):
    """Return the labels and values of the replicate whose cross has the scale of ``frames[0]``."""
    atoms.universe.trajectory[int(frames[0])]
    scale = round(float(atoms.positions[0, 0]), 3)
    CALLS.append(scale)
    labels, values = TABLE[scale]
    return list(labels), np.asarray(values, dtype=float)


def _wrong_shape(atoms, frames):
    return ["a", "b", "c"], np.zeros(2)


@pytest.fixture(autouse=True)
def _reset_calls():
    CALLS.clear()
    yield
    CALLS.clear()


def _config(tmp_path: Path, label: str, scales) -> Path:
    config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
    for replicate, scale in enumerate(scales, start=1):
        write_openmm_replicate(config, replicate, [scale] * 4)
    return config


@pytest.fixture()
def study(tmp_path: Path):
    configs = {label: _config(tmp_path, label, scales) for label, scales in SCALES.items()}
    return pz.Study.from_configs(configs, equilibration="0ns")


def _run(study, tmp_path, **options):
    return study.per_replicate(
        _measure,
        pz.select("all"),
        unit=None,
        labels="returned",
        name="pairs",
        output_dir=tmp_path / "out",
        **{"missing": 0.0, **options},
    )


def _as_dicts(values) -> dict[str, list[dict]]:
    return {
        condition: [dict(zip(values.labels, row.tolist())) for row in rows]
        for condition, rows in values.values.items()
    }


def _folder(tmp_path: Path, condition: str, replicate: int) -> Path:
    return tmp_path / "out" / "polyzymd_results" / "pairs" / condition / f"replicate_{replicate}"


def _expected(condition: str) -> list[dict]:
    order = ["12-SBM", "45-EGM", "7-SBM"]
    return [
        {label: dict(zip(*TABLE[scale])).get(label, 0.0) for label in order}
        for scale in SCALES[condition]
    ]


def test_replicates_with_different_labels_line_up_with_missing_filled(study, tmp_path) -> None:
    values = _run(study, tmp_path)

    # Labels in the order they first appear, across every replicate of every condition.
    assert values.labels == ["12-SBM", "45-EGM", "7-SBM"]
    assert _as_dicts(values) == {"A": _expected("A"), "B": _expected("B")}
    assert sorted(CALLS) == [1.0, 1.0, 2.0, 3.0, 4.0, 5.0]


def test_the_missing_value_fills_every_absent_label(study, tmp_path) -> None:
    values = _run(study, tmp_path, missing=-1.0)

    replicate_3 = dict(zip(values.labels, values.values["A"][2].tolist()))
    assert replicate_3 == {"12-SBM": -1.0, "45-EGM": -1.0, "7-SBM": -1.0}


def test_each_replicate_stores_its_own_labels_beside_its_values(study, tmp_path) -> None:
    import json

    _run(study, tmp_path)

    for condition, scales in SCALES.items():
        for replicate, scale in enumerate(scales, start=1):
            folder = _folder(tmp_path, condition, replicate)
            assert json.loads((folder / "labels.json").read_text()) == TABLE[scale][0]
            with np.load(folder / "values.npz") as data:
                assert data["values"].tolist() == pytest.approx(TABLE[scale][1])
            assert json.loads((folder / "record.json").read_text())["labels"] == "returned"


def test_stored_results_are_reused_without_calling_the_function(study, tmp_path) -> None:
    first = _run(study, tmp_path)
    stored = sorted((tmp_path / "out" / "polyzymd_results" / "pairs").rglob("*.*"))
    before = [path.stat().st_mtime_ns for path in stored]
    CALLS.clear()

    again = _run(study, tmp_path)

    assert CALLS == []
    assert [path.stat().st_mtime_ns for path in stored] == before
    assert again.labels == first.labels
    assert _as_dicts(again) == _as_dicts(first)


def test_a_missing_labels_file_forces_that_replicate_to_be_recomputed(study, tmp_path) -> None:
    first = _run(study, tmp_path)
    (_folder(tmp_path, "A", 2) / "labels.json").unlink()
    CALLS.clear()

    again = _run(study, tmp_path)

    assert CALLS == [2.0]
    assert (_folder(tmp_path, "A", 2) / "labels.json").is_file()
    assert _as_dicts(again) == _as_dicts(first)


def test_a_corrupt_labels_file_forces_that_replicate_to_be_recomputed(study, tmp_path) -> None:
    _run(study, tmp_path)
    (_folder(tmp_path, "B", 1) / "labels.json").write_text("{not json")
    CALLS.clear()

    again = _run(study, tmp_path)

    assert CALLS == [4.0]
    assert _as_dicts(again)["B"] == _expected("B")


def test_recompute_calls_the_function_for_every_replicate(study, tmp_path) -> None:
    _run(study, tmp_path)
    CALLS.clear()

    _run(study, tmp_path, recompute=True)

    assert len(CALLS) == 6


def test_labels_given_as_a_list_do_not_reuse_returned_labels(study, tmp_path) -> None:
    """The record holds ``"returned"``, so a call with fixed labels does not match it."""
    _run(study, tmp_path)
    CALLS.clear()

    def fixed(atoms, frames):
        return np.zeros(3)

    study.per_replicate(
        fixed,
        pz.select("all"),
        unit=None,
        labels=["12-SBM", "45-EGM", "7-SBM"],
        name="pairs",
        output_dir=tmp_path / "out",
    )
    again = _run(study, tmp_path)

    assert sorted(CALLS) == [1.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    assert _as_dicts(again) == {"A": _expected("A"), "B": _expected("B")}


def test_without_missing_different_labels_are_refused(study, tmp_path) -> None:
    with pytest.raises(ProtocolError, match="has no value for labels") as err:
        _run(study, tmp_path, missing=None)

    assert "missing=" in err.value.hint


def test_a_function_returning_values_that_do_not_match_its_labels_is_refused(
    study, tmp_path
) -> None:
    with pytest.raises(ProtocolError, match=r"returned shape \(2,\) for 3 labels"):
        study.per_replicate(
            _wrong_shape,
            pz.select("all"),
            unit=None,
            labels="returned",
            missing=0.0,
            output_dir=tmp_path / "out",
        )


def test_numpy_labels_are_stored_as_plain_json_numbers(tmp_path) -> None:
    import json

    config = _config(tmp_path, "C", (6.0,))
    study = pz.Study.from_configs({"C": config}, equilibration="0ns")

    values = study.per_replicate(
        _measure, pz.select("all"), unit=None, labels="returned", output_dir=tmp_path / "out"
    )
    CALLS.clear()
    again = study.per_replicate(
        _measure, pz.select("all"), unit=None, labels="returned", output_dir=tmp_path / "out"
    )

    folder = tmp_path / "out" / "polyzymd_results" / "measure" / "C" / "replicate_1"
    assert json.loads((folder / "labels.json").read_text()) == [3, 1]
    assert values.labels == again.labels == [3, 1]
    assert CALLS == []
    assert values.values["C"][0].tolist() == pytest.approx([0.4, 0.6])


def test_returned_labels_compare_label_by_label(study, tmp_path) -> None:
    values = _run(study, tmp_path)

    report = values.compare()

    entries = {row.entry for row in report.pairwise}
    assert entries == {"12-SBM", "45-EGM", "7-SBM"}
