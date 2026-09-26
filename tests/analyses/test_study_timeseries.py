"""Known-answer tests for Study, Study.timeseries and the replicate statistics.

Every replicate is a real OpenMM run directory with one DCD segment of four
unit-mass atoms on a cross. Frame ``k`` is scaled by ``scales[k]``, so its
radius of gyration is ``scales[k]`` and its time is ``k * 0.1`` ns.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg
from polyzymd.analyses.shared.statistics import mean_sem_ci
from tests._support.analysis_testkit import (
    replicate_values,
    write_openmm_replicate,
    write_simulation_config,
)

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

# Frame k is at 0.1 * k ns, and a 0.25 ns window leaves frames 3 to 9.
EQUILIBRATION = "0.25ns"
FRAMES = list(range(3, 10))


def radius_of_gyration(atoms):
    """Mass-weighted radius of gyration of ``atoms``."""
    return atoms.radius_of_gyration()


def _scales(base: float) -> list[float]:
    """Ten per-frame Rg values whose production mean is ``base + 0.06``."""
    return [base + 0.01 * k for k in range(10)]


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    """Two conditions of three replicates, B larger than A by one Å."""
    paths = {}
    for label, offset in (("A", 1.0), ("B", 2.0)):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            write_openmm_replicate(config, replicate, _scales(offset + 0.1 * replicate))
        paths[label] = config
    return paths


@pytest.fixture()
def study(configs: dict[str, Path]) -> pz.Study:
    """The two-condition study with the 0.25 ns window."""
    return pz.Study.from_configs(configs, equilibration=EQUILIBRATION)


class TestStudy:
    """Conditions, replicates, frames, times and identity."""

    def test_labels_control_and_iteration(self, study: pz.Study) -> None:
        assert study.labels == ["A", "B"]
        assert study.control == "A"
        assert [condition.label for condition in study] == ["A", "B"]
        assert [r.index for r in study["B"].replicates] == [1, 2, 3]

    def test_sequence_labels_come_from_folders(self, configs: dict[str, Path]) -> None:
        assert pz.Study.from_configs(list(configs.values()), equilibration="0ns").labels == [
            "A",
            "B",
        ]

    def test_frames_and_times_follow_the_window(self, study: pz.Study) -> None:
        replicate = study["A"].replicates[0]
        assert replicate.frames.tolist() == FRAMES
        assert replicate.times == pytest.approx([0.1 * k for k in FRAMES], abs=1e-6)

    def test_universe_is_cached(self, study: pz.Study) -> None:
        replicate = study["A"].replicates[0]
        assert replicate.universe() is replicate.universe()

    def test_identity_records_config_window_and_files(self, study: pz.Study) -> None:
        identity = study["A"].replicates[0].identity
        assert identity["config_hash"] == study["A"].config_hash
        assert identity["equilibration"] == EQUILIBRATION
        assert identity["topology"]["path"].endswith("solvated_system.pdb")
        assert identity["trajectories"][0]["path"].endswith("production_0_trajectory.dcd")
        assert identity["trajectories"][0]["size_bytes"] > 0

    def test_replicate_subset_and_missing_replicate(self, configs: dict[str, Path]) -> None:
        subset = pz.Study.from_configs(configs, equilibration="0ns", replicates=[2])
        assert [r.index for r in subset["A"].replicates] == [2]
        with pytest.raises(ProtocolError, match="no run directory"):
            pz.Study.from_configs(configs, equilibration="0ns", replicates=[7])

    def test_unknown_label_and_bad_window(self, study: pz.Study, configs) -> None:
        with pytest.raises(ProtocolError, match="no condition"):
            study["C"]
        with pytest.raises(ProtocolError, match="equilibration"):
            pz.Study.from_configs(configs, equilibration="soon")

    def test_importing_polyzymd_imports_no_numpy(self) -> None:
        import subprocess

        code = "import sys, polyzymd; print('numpy' in sys.modules or 'MDAnalysis' in sys.modules)"
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        assert result.stdout.strip() == "False"


def _run(study: pz.Study, tmp_path: Path, function=radius_of_gyration, selection="all", **kw):
    """Run the Rg timeseries into ``tmp_path``."""
    return study.timeseries(
        function, pz.select(selection), unit="A", name="rg", output_dir=tmp_path, **kw
    )


class TestTimeseries:
    """Per-frame values, storage and reuse."""

    def test_values_are_the_known_radii(self, study: pz.Study, tmp_path: Path) -> None:
        series = _run(study, tmp_path).series["A"][0]
        assert series.values == pytest.approx([1.1 + 0.01 * k for k in FRAMES], abs=1e-5)
        assert series.frames.tolist() == FRAMES

    def test_record_holds_code_arguments_inputs_and_versions(
        self, study: pz.Study, tmp_path: Path
    ) -> None:
        folder = tmp_path / "polyzymd_results" / "rg" / "A" / "replicate_1"
        _run(study, tmp_path)
        record = json.loads((folder / "record.json").read_text())
        assert record["function"]["qualname"] == "radius_of_gyration"
        assert record["function"]["module"] == __name__
        assert record["function"]["hash_of"] == "source"
        assert record["arguments"]["args"] == [{"select": "all"}]
        assert record["config_hash"] == study["A"].config_hash
        assert record["equilibration"] == EQUILIBRATION
        assert record["frames"] == FRAMES
        assert len(record["times_ns"]) == len(FRAMES)
        assert record["trajectories"][0]["path"].endswith(".dcd")
        assert record["unit"] == "A"
        assert set(record["versions"]) == {"polyzymd", "MDAnalysis", "numpy", "python"}
        with np.load(folder / "series.npz") as data:
            assert data["frames"].tolist() == FRAMES

    def test_matching_record_is_reused(self, study: pz.Study, tmp_path, monkeypatch) -> None:
        _run(study, tmp_path)
        monkeypatch.setattr(
            "MDAnalysis.analysis.base.AnalysisFromFunction.run",
            lambda *args, **kwargs: pytest.fail("a matching stored series was measured again"),
        )
        assert _run(study, tmp_path).series["B"][2].values[0] == pytest.approx(2.33, abs=1e-5)

    def test_recompute_measures_again(self, study: pz.Study, tmp_path: Path) -> None:
        folder = tmp_path / "polyzymd_results" / "rg" / "A" / "replicate_1"
        _run(study, tmp_path)
        (folder / "series.npz").unlink()
        np.savez(folder / "series.npz", values=np.zeros(7))
        assert _run(study, tmp_path).series["A"][0].values[0] == 0.0
        assert _run(study, tmp_path, recompute=True).series["A"][0].values[0] > 1.0

    def test_changed_argument_recomputes(self, study: pz.Study, tmp_path: Path) -> None:
        _run(study, tmp_path)
        # Two opposite atoms of the cross have the same radius of gyration.
        values = _run(study, tmp_path, selection="name C1 C2").series["A"][0].values
        record = json.loads(
            (tmp_path / "polyzymd_results/rg/A/replicate_1/record.json").read_text()
        )
        assert record["arguments"]["args"] == [{"select": "name C1 C2"}]
        assert values == pytest.approx([1.1 + 0.01 * k for k in FRAMES], abs=1e-5)

    def test_changed_function_source_recomputes(
        self, study: pz.Study, tmp_path: Path, monkeypatch
    ) -> None:
        module_dir = tmp_path / "code"
        module_dir.mkdir()
        source = module_dir / "measure_rg.py"
        source.write_text("def measure(atoms):\n    return atoms.radius_of_gyration()\n")
        monkeypatch.syspath_prepend(str(module_dir))
        module = importlib.import_module("measure_rg")
        first = _run(study, tmp_path, function=module.measure).series["A"][0].values
        source.write_text("def measure(atoms):\n    return 2 * atoms.radius_of_gyration()\n")
        module = importlib.reload(module)
        second = _run(study, tmp_path, function=module.measure).series["A"][0].values
        sys.modules.pop("measure_rg", None)
        assert second == pytest.approx(2 * first)

    def test_changed_input_file_recomputes(self, configs, tmp_path: Path) -> None:
        study = pz.Study.from_configs(configs, equilibration=EQUILIBRATION)
        _run(study, tmp_path)
        write_openmm_replicate(configs["A"], 1, [5.0] * 10)
        fresh = pz.Study.from_configs(configs, equilibration=EQUILIBRATION)
        assert _run(fresh, tmp_path).series["A"][0].values == pytest.approx([5.0] * 7)

    def test_array_output_is_rejected(self, study: pz.Study, tmp_path: Path) -> None:
        with pytest.raises(ProtocolError, match="one number per frame"):
            study.timeseries(lambda atoms: atoms.positions[0], pz.select("all"), unit=None)

    def test_universe_argument(self, study: pz.Study, tmp_path: Path) -> None:
        def n_atoms(u):
            return len(u.atoms)

        series = study.timeseries(n_atoms, pz.universe(), unit=None, output_dir=tmp_path)
        assert series.series["A"][0].values.tolist() == [4.0] * 7


class TestReduce:
    """Named reductions and a function reduction."""

    def test_named_and_callable_reductions(self, study: pz.Study, tmp_path: Path) -> None:
        series = _run(study, tmp_path)
        assert series.reduce("mean").values["A"] == pytest.approx([1.16, 1.26, 1.36], abs=1e-5)
        expected_std = float(np.std([0.01 * k for k in FRAMES], ddof=1))
        assert series.reduce("std").values["B"] == pytest.approx([expected_std] * 3, abs=1e-5)
        last = series.reduce(lambda values, times: float(values[-1]))
        assert last.values["A"] == pytest.approx([1.19, 1.29, 1.39], abs=1e-5)

    def test_fraction_needs_zeros_and_ones(self, study: pz.Study, tmp_path: Path) -> None:
        with pytest.raises(ProtocolError, match="0 and 1"):
            _run(study, tmp_path).reduce("fraction")


class TestSummaryAndCompare:
    """Intervals, tests and warnings with known answers."""

    def test_summary_matches_mean_sem_ci(self) -> None:
        report = replicate_values({"A": [1.0, 2.0, 4.0]}).summary()
        (row,) = report.conditions
        expected = mean_sem_ci([1.0, 2.0, 4.0])
        assert (row.mean, row.sem, row.ci95) == (
            expected.mean,
            expected.sem,
            (expected.ci_low, expected.ci_high),
        )
        assert row.replicate_values == [1.0, 2.0, 4.0]
        assert row.statistical_inefficiency == [1.0, 1.0, 1.0]
        assert row.n_effective == [20.0, 20.0, 20.0]

    def test_identical_replicates_are_not_estimable(self) -> None:
        report = replicate_values({"A": [3.0, 3.0, 3.0]}).summary()
        assert report.conditions[0].ci95 is None
        assert report.conditions[0].ci_method == "not_estimable"
        assert any("not estimable" in text for text in report.warnings)

    def test_fraction_interval_past_one_is_flagged(self) -> None:
        report = replicate_values({"A": [1.0, 1.0, 0.0]}, how="fraction").summary()
        assert report.unit is None
        assert any("fraction bounds" in text for text in report.warnings)

    def test_compare_welch_student_bh_and_cohens_d(self) -> None:
        data = {"A": [1.0, 1.2, 1.1], "B": [2.0, 2.4, 2.2], "C": [1.0, 1.3, 1.15]}
        values = replicate_values(data)
        welch = values.compare()
        student = values.compare(test="student")
        for row in welch.pairwise:
            expected = stats.ttest_ind(data[row.b], data[row.a], equal_var=False)
            assert row.p == pytest.approx(expected.pvalue)
            low, high = expected.confidence_interval(0.95)
            assert row.delta_ci95 == pytest.approx((low, high))
        assert student.pairwise[0].p == pytest.approx(stats.ttest_ind(data["B"], data["A"]).pvalue)
        adjusted = benjamini_hochberg([row.p for row in welch.pairwise])
        assert [row.p_adjusted for row in welch.pairwise] == [
            item.adjusted_p_value for item in adjusted
        ]
        # B is larger than A, so the difference and Cohen's d are both positive.
        assert welch.pairwise[0].delta > 0 and welch.pairwise[0].cohens_d > 0
        assert welch.pairwise[0].test == "welch_t"

    def test_constant_conditions_are_not_testable(self) -> None:
        report = replicate_values({"A": [0.0, 0.0], "B": [0.0, 0.0], "C": [1.0, 2.0]}).compare()
        untestable, testable = report.pairwise
        assert not untestable.testable and untestable.p is None
        assert untestable.p_adjusted is None
        # The untestable row takes no part in the correction.
        assert testable.p_adjusted == pytest.approx(testable.p)

    def test_conditions_and_control_arguments(self) -> None:
        values = replicate_values({"A": [1.0, 2.0], "B": [3.0, 4.0], "C": [5.0, 6.0]})
        report = values.compare(control="B", conditions=["C"])
        assert [(row.a, row.b) for row in report.pairwise] == [("B", "C")]
        assert [row.label for row in values.summary(conditions=["C"]).conditions] == ["C"]

    def test_agent_text_lists_every_value(self) -> None:
        text = (
            replicate_values({"A": [1.0, 2.0, 4.0], "B": [2.0, 3.0, 5.0]}).compare().to_agent_text()
        )
        assert "values 1, 2, 4" in text and "g 1, 1, 1" in text and "A vs B" in text


def test_statistics_match_the_stored_legacy_rg_comparison() -> None:
    """Replicate means of a real legacy rg comparison give its SEM, intervals and p values.

    ``tests/data/comparison_artifacts/rg.json`` is the comparison the legacy rg
    plugin stored for a LipA campaign, with Student t tests and one
    Benjamini-Hochberg family over the 15 comparisons of its two runs.
    """
    from tests.analyses.test_protocols_real_artifacts import _report

    stored = json.loads(
        (Path(__file__).parent.parent / "data/comparison_artifacts/rg.json").read_text()
    )
    raw: dict[tuple[str, str, str], float] = {}
    for run in stored["run_labels"]:
        per_condition = {
            condition["label"]: summary["per_replicate_means"]
            for condition in stored["conditions"]
            for summary in condition["run_summaries"]
            if summary["label"] == run
        }
        sems = {
            condition["label"]: summary["sem_rg"]
            for condition in stored["conditions"]
            for summary in condition["run_summaries"]
            if summary["label"] == run
        }
        legacy = {row.label: row for row in _report("rg", run=run).conditions}
        values = replicate_values(per_condition)
        for row in values.summary().conditions:
            assert row.sem == pytest.approx(sems[row.label], rel=1e-12)
            assert row.mean == pytest.approx(legacy[row.label].mean, rel=1e-12)
            assert row.ci95 == pytest.approx(legacy[row.label].ci95, rel=1e-12)
        labels = list(per_condition)
        for position, control in enumerate(labels):
            report = values.compare(control=control, conditions=labels[position:], test="student")
            raw.update({(run, row.a, row.b): row.p for row in report.pairwise})
    rows = stored["pairwise_comparisons"]
    family = benjamini_hochberg(
        [raw[(r["run_label"], r["condition_a"], r["condition_b"])] for r in rows]
    )
    for row, corrected in zip(rows, family, strict=True):
        key = (row["run_label"], row["condition_a"], row["condition_b"])
        assert raw[key] == pytest.approx(row["p_value"], rel=1e-9)
        assert corrected.adjusted_p_value == pytest.approx(row["p_value_adjusted"], rel=1e-9)
