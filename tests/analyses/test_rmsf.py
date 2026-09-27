"""Tests for labelled per-replicate results, the rmsf function and ``polyzymd analyze rmsf``.

Each replicate is an OpenMM run directory with one DCD segment of three
residues of two atoms each. Every frame is the same shape with atom noise
whose size grows with the residue, rotated and translated at random, so the
fluctuation is known to be largest at residue 3. Frame ``k`` is at
``0.1 * k`` ns, and the 0.25 ns window leaves frames 3 to 9.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import rmsf
from polyzymd.analyses.shared import plotting
from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg
from polyzymd.analyses.shared.statistics import mean_sem_ci
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0.25ns"
RESINDEX = [0, 0, 1, 1, 2, 2]
SHAPE = np.array(
    [[0, 0, 0], [1.5, 0, 0], [3, 1, 0], [4.5, 1, 0], [6, 0, 1], [7.5, 0, 1]], dtype=float
)


def _frames(seed: int, noise: float, n_frames: int = 10) -> np.ndarray:
    """Noisy, rotated and translated copies of SHAPE, noise growing with the residue."""
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(seed)
    scale = noise * np.repeat([0.2, 0.5, 1.0], 2)[:, np.newaxis]
    return np.array(
        [
            Rotation.random(random_state=rng).apply(SHAPE + scale * rng.normal(size=SHAPE.shape))
            + rng.normal(scale=5.0, size=3)
            for _ in range(n_frames)
        ]
    )


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    """Two conditions of three replicates, B fluctuating twice as much as A."""
    paths = {}
    for label, noise in (("A", 0.3), ("B", 0.6)):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            write_openmm_frames(
                config, replicate, _frames(10 * replicate + len(label), noise), RESINDEX
            )
        paths[label] = config
    return paths


@pytest.fixture()
def study(configs) -> pz.Study:
    return pz.Study.from_configs(configs, equilibration=EQUILIBRATION)


def _profile(study, tmp_path, mode="average", **options):
    return study.per_replicate(
        rmsf,
        pz.select("all"),
        pz.select("all"),
        pz.reference(mode, "(all) or (all)", alignment="all", **options),
        unit="A",
        labels=lambda u: u.residues.resids,
        output_dir=tmp_path,
    )


def _by_hand(coordinates: np.ndarray, reference: np.ndarray) -> list[np.ndarray]:
    """Align with AlignTraj onto ``reference``; return per-atom deviation, rms.RMSF and offset."""
    from MDAnalysis.analysis import align, rms

    mobile = mda.Universe.empty(6, n_residues=3, atom_resindex=RESINDEX, trajectory=True)
    mobile.load_new(coordinates.astype(np.float32), format="MEMORY")
    ref = mda.Merge(mobile.atoms)
    ref.load_new(reference.astype(np.float32)[np.newaxis], format="MEMORY")
    align.AlignTraj(mobile, ref, select="all", in_memory=True).run()
    moved = np.array([mobile.atoms.positions for _ in mobile.trajectory], float)
    deviation = np.sqrt(np.mean(np.sum((moved - reference) ** 2, axis=2), axis=0))
    offset = np.linalg.norm(moved.mean(axis=0) - reference, axis=1)
    return [deviation, rms.RMSF(mobile.atoms).run().results.rmsf, offset]


def _residues(per_atom: np.ndarray) -> np.ndarray:
    return per_atom.reshape(3, 2).mean(axis=1)


def test_rmsf_deviation_and_offset_equal_aligntraj_and_gmx_definitions() -> None:
    """rmsf is rms.RMSF after AlignTraj, rms_deviation the gmx rmsf -od deviation, in one pass."""
    from polyzymd.analyses.functions import (
        RMS_PARTS,
        _superposed_deviations,
        rms_decomposition,
        rms_deviation,
    )

    coordinates = _frames(3, 0.5)
    universe = mda.Universe.empty(6, n_residues=3, atom_resindex=RESINDEX, trajectory=True)
    universe.load_new(coordinates.astype(np.float32), format="MEMORY")
    before = universe.trajectory.coordinate_array.copy()
    ref = mda.Merge(universe.atoms)
    ref.load_new(coordinates[4][np.newaxis].astype(np.float32), format="MEMORY")
    frames = np.arange(2, 10)
    expected = _by_hand(coordinates[2:], coordinates[4].astype(np.float32).astype(float))
    atoms = universe.atoms
    assert rmsf(atoms, atoms, ref.atoms, frames) == pytest.approx(_residues(expected[1]), abs=1e-5)
    assert rms_deviation(atoms, atoms, ref.atoms, frames) == pytest.approx(
        _residues(expected[0]), abs=1e-5
    )
    assert RMS_PARTS == ("rms_deviation", "rmsf", "offset")
    both = rms_decomposition(atoms, atoms, ref.atoms, frames)
    assert both[:3] == pytest.approx(np.vstack([_residues(e) for e in expected]), abs=1e-5)
    assert both[3:] == pytest.approx(np.vstack([_residues(e**2) for e in expected]), abs=1e-5)
    assert np.max(np.abs(both[3] - both[4] - both[5])) < 1e-5
    deviation, fluctuation, offset = _superposed_deviations(atoms, atoms, ref.atoms, frames)
    assert np.max(np.abs(deviation**2 - fluctuation**2 - offset**2)) < 1e-5
    # Fitting on residues 1 and 2 only gives other values than fitting on all atoms.
    assert not np.allclose(rmsf(atoms, atoms[:4], ref.atoms, frames), both[1], atol=1e-3)
    assert np.array_equal(universe.trajectory.coordinate_array, before)
    with pytest.raises(ProtocolError, match="reference has 4 atoms"):
        rmsf(atoms, atoms, ref.atoms[:4], frames)


def test_per_replicate_profile_is_labelled_by_residue(study, tmp_path) -> None:
    """Each replicate's profile equals the hand-aligned RMSF, and summary rows are per residue."""
    profile = _profile(study, tmp_path, "frame", frame=1)
    assert profile.labels == [1, 2, 3]
    replicate = study["B"].replicates[1]
    u = replicate.universe()
    coordinates = np.array([u.atoms.positions for _ in u.trajectory[replicate.frames]], float)
    expected = _residues(_by_hand(coordinates, coordinates[0])[1])
    assert profile.values["B"][1] == pytest.approx(expected, abs=1e-5)
    assert expected[2] > expected[0]
    summary = profile.summary()
    assert [(row.label, row.entry) for row in summary.conditions] == [
        ("A", "1"),
        ("B", "1"),
        ("A", "2"),
        ("B", "2"),
        ("A", "3"),
        ("B", "3"),
    ]
    row = summary.conditions[5]
    column = [values[2] for values in profile.values["B"]]
    stats = mean_sem_ci(column)
    assert row.replicate_values == pytest.approx(column)
    assert (row.mean, row.sem, *row.ci95) == pytest.approx(
        (stats.mean, stats.sem, stats.ci_low, stats.ci_high)
    )
    text = summary.to_agent_text().split("\n")
    assert (
        text[0].startswith("# polyzymd analyze rmsf  metric rmsf  unit A")
        and "conditions 2" in text[0]
    )
    assert text[6].startswith("3  B  n 3  mean ")


def test_profile_record_holds_the_labels_and_the_reference_hash(study, tmp_path) -> None:
    import json

    reference = tmp_path / "ref.pdb"
    universe = study["A"].replicates[0].universe()
    universe.trajectory[0]
    universe.atoms.write(str(reference))
    profile = _profile(study, tmp_path, "external", file=reference)
    record = json.loads((profile.source.series["A"][0] / "record.json").read_text())
    assert record["labels"] == [1, 2, 3]
    assert len(record["arguments"]["args"][2]["reference"]["file"]["sha256"]) == 64
    again = _profile(study, tmp_path, "external", file=reference)
    assert again.values["A"][0] == pytest.approx(profile.values["A"][0])


def test_mean_rmsf_is_the_mean_over_residues(study, tmp_path) -> None:
    profile = _profile(study, tmp_path)
    mean = profile.over_labels("mean")
    assert mean.metric == "mean_rmsf" and mean.labels is None
    for label in ("A", "B"):
        assert mean.values[label] == pytest.approx(
            [float(np.mean(v)) for v in profile.values[label]]
        )
    report = mean.compare()
    assert [row.entry for row in report.pairwise] == [None]
    assert report.verdict[0].startswith("B larger mean_rmsf than A")
    with pytest.raises(ProtocolError, match="labelled values"):
        mean.over_labels()


def test_per_label_compare_corrects_over_every_label_and_condition(study, tmp_path) -> None:
    """One Benjamini-Hochberg family covers every residue of every compared condition."""
    profile = _profile(study, tmp_path)
    report = profile.compare()
    assert [row.entry for row in report.pairwise] == ["1", "2", "3"]
    assert {row.family_size for row in report.pairwise} == {3}
    expected = benjamini_hochberg([row.p for row in report.pairwise])
    assert [row.p_adjusted for row in report.pairwise] == pytest.approx(
        [item.adjusted_p_value for item in expected]
    )
    assert report.pairwise[0].test == "welch_t"
    assert report.verdict[0].startswith("B vs A: ")
    skipped = profile.compare(untested=[2])
    assert [row.entry for row in skipped.pairwise] == ["1", "3"]
    assert {row.family_size for row in skipped.pairwise} == {2}
    assert any("left out of the tests" in text for text in skipped.warnings)
    assert len(skipped.conditions) == 6


def _resid_value(u, frames):
    """Ten times each residue ID, whatever the order of the residues."""
    return u.residues.resids * 10.0


def test_labels_align_across_replicates_and_missing_labels(tmp_path) -> None:
    """Replicates are lined up by label; a missing label fails unless missing= fills it."""
    config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "A" / "scratch")
    write_openmm_frames(config, 1, _frames(1, 0.3), RESINDEX)
    write_openmm_frames(config, 2, _frames(2, 0.3), RESINDEX, resids=[3, 2, 1])
    study = pz.Study.from_configs({"A": config}, equilibration="0ns")
    labels = lambda u: u.residues.resids  # noqa: E731
    values = study.per_replicate(
        _resid_value, pz.universe(), unit=None, labels=labels, output_dir=tmp_path
    )
    assert values.labels == [1, 2, 3]
    assert [list(v) for v in values.values["A"]] == [[10.0, 20.0, 30.0]] * 2
    write_openmm_frames(config, 3, _frames(3, 0.3)[:, :4], RESINDEX[:4])
    study = pz.Study.from_configs({"A": config}, equilibration="0ns")
    with pytest.raises(ProtocolError, match=r"replicate 3 has no value for labels \[3\]"):
        study.per_replicate(
            _resid_value, pz.universe(), unit=None, labels=labels, output_dir=tmp_path
        )
    filled = study.per_replicate(
        _resid_value, pz.universe(), unit=None, labels=labels, missing=0.0, output_dir=tmp_path
    )
    assert list(filled.values["A"][2]) == [10.0, 20.0, 0.0]
    with pytest.raises(ProtocolError, match="returned shape"):
        study.per_replicate(
            _resid_value,
            pz.universe(),
            unit=None,
            labels=[1, 2],
            output_dir=tmp_path,
            recompute=True,
        )


def test_constant_labels_are_not_testable_and_leave_the_family(tmp_path) -> None:
    config_a = write_simulation_config(tmp_path / "A", scratch=tmp_path / "A" / "scratch")
    config_b = write_simulation_config(tmp_path / "B", scratch=tmp_path / "B" / "scratch")
    for config in (config_a, config_b):
        for replicate in (1, 2):
            write_openmm_frames(config, replicate, _frames(replicate, 0.3), RESINDEX)
    study = pz.Study.from_configs({"A": config_a, "B": config_b}, equilibration="0ns")
    values = study.per_replicate(
        _resid_value,
        pz.universe(),
        unit=None,
        labels=lambda u: u.residues.resids,
        output_dir=tmp_path,
    )
    report = values.compare()
    assert all(not row.testable and row.family_size is None for row in report.pairwise)
    assert sum("is not testable" in text for text in report.warnings) == 3


@pytest.fixture()
def figures(monkeypatch) -> dict[str, object]:
    """Keep every saved figure open, by file stem, so its Axes can be read."""
    saved: dict[str, object] = {}
    original = plotting.save_figure

    def keep(fig, output_path, plot_settings, **kwargs):
        saved[Path(output_path).stem] = fig
        return original(fig, output_path, plot_settings, close=False)

    monkeypatch.setattr(plotting, "save_figure", keep)
    return saved


def test_cli_rmsf_reports_mean_rmsf_and_draws_the_profile(configs, tmp_path, figures) -> None:
    """The profile lines and band equal the stored per-residue values and intervals."""
    arguments = ["rmsf", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--set", "selection=all", "--set", "alignment_selection=all"]
    arguments += ["--set", "reference_mode=average", "--set", "highlight_residues=[2]"]
    arguments += ["--output-dir", str(tmp_path)]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    lines = result.stdout.strip().split("\n")
    assert lines[0].startswith("# polyzymd analyze rmsf  metric core_rmsf  unit A  run core_rmsf")
    assert lines[-1].startswith("verdict: B larger core_rmsf than A")
    folder = tmp_path / "figures" / "rmsf"
    assert {path.name for path in folder.iterdir()} == {
        "rmsf_profile.png",
        "rms_deviation_profile.png",
        "offset_profile.png",
        "rms_decomposition.png",
        "rmsf_comparison.png",
    }
    decomposition = figures["rms_decomposition"].axes
    assert [ax.get_title() for ax in decomposition] == ["A (n = 3)", "B (n = 3)"]
    assert [line.get_label() for line in decomposition[0].get_lines()] == [
        "rms deviation",
        "rmsf",
        "offset",
    ]
    per_residue = CliRunner().invoke(analyze_command, [*arguments, "--run", "rmsf"])
    assert per_residue.exit_code == 0, per_residue.output
    rows = [line for line in per_residue.stdout.split("\n") if "  A vs B  " in line]
    assert [row.split()[0] for row in rows] == ["1", "2", "3"] and "family 3" in rows[0]
    stored = _profile(pz.Study.from_configs(configs, equilibration=EQUILIBRATION), tmp_path)
    summary = stored.summary()
    ax = figures["rmsf_profile"].axes[0]
    for label in ("A", "B"):
        rows = [row for row in summary.conditions if row.label == label]
        mean = next(line for line in ax.get_lines() if line.get_label() == f"{label} (n = 3)")
        assert list(mean.get_xdata()) == [1, 2, 3]
        assert mean.get_ydata() == pytest.approx([row.mean for row in rows])
        for values in stored.values[label]:
            assert any(np.allclose(line.get_ydata(), values) for line in ax.get_lines())
    bands = [item.get_paths()[0].vertices for item in ax.collections]
    for label, band in zip(("A", "B"), bands, strict=True):
        for residue, row in zip((1, 2, 3), [r for r in summary.conditions if r.label == label]):
            edge = sorted(band[np.isclose(band[:, 0], residue), 1])
            assert (edge[0], edge[-1]) == pytest.approx(tuple(row.ci95))
    highlight = [line for line in ax.get_lines() if line.get_color() == "red"]
    assert list(highlight[0].get_xdata()) == [2, 2]
    none = tmp_path / "none"
    quiet = CliRunner().invoke(analyze_command, [*arguments[:-1], str(none), "--no-plots"])
    assert quiet.exit_code == 0 and not (none / "figures").exists()


def test_analyze_rmsf_refuses_an_unknown_run(configs) -> None:
    from polyzymd.analyses import analyze

    with pytest.raises(ProtocolError, match="no result named 'x'"):
        analyze(
            "rmsf",
            [configs["A"]],
            equilibration=EQUILIBRATION,
            run="x",
            settings={"selection": "all", "alignment_selection": "all"},
        )


def test_rms_deviation_defaults_to_the_reference_file_and_reports_every_part(
    configs, study, tmp_path
) -> None:
    """With a reference file, rms_deviation superposes on it; runs pick any of the six results."""
    from polyzymd.analyses import analyze

    reference = tmp_path / "ref.pdb"
    universe = study["A"].replicates[0].universe()
    universe.trajectory[0]
    universe.atoms.write(str(reference))
    settings = {"selection": "all", "alignment_selection": "all", "reference_file": str(reference)}
    options = {"equilibration": EQUILIBRATION, "settings": settings, "output_dir": tmp_path}
    report = analyze("rms_deviation", [configs["A"], configs["B"]], plots=False, **options)
    assert report.analysis == "rms_deviation" and report.metric == "core_rms_deviation"
    parts = ("rms_deviation", "rmsf", "offset")
    assert report.all_runs == [f"{kind}_{p}" for kind in ("core", "mean") for p in parts] + list(
        parts
    )
    rows = study.per_replicate(
        functions.rms_decomposition,
        pz.select("all"),
        pz.select("all"),
        pz.reference("external", "(all) or (all)", file=reference, alignment="all"),
        unit="A",
        labels=lambda u: u.residues.resids,
        output_dir=tmp_path / "check",
        parts=functions.RMS_PARTS + functions.MS_PARTS,
    )
    core = [float(np.sqrt(np.mean(v))) for v in rows["ms_deviation"].values["A"]]
    assert report.conditions[0].replicate_values == pytest.approx(core)
    assert report.provenance.settings["residues"] == {"core": [1, 2, 3]}
    assert report.provenance.settings["reference_mode"] == "external"
    offset = analyze("rms_deviation", [configs["A"]], plots=False, run="offset", **options)
    assert [row.entry for row in offset.conditions] == ["1", "2", "3"]
    mean = analyze("rms_deviation", [configs["A"]], plots=False, run="mean_rmsf", **options)
    assert mean.conditions[0].replicate_values == pytest.approx(
        [float(np.mean(v)) for v in rows["rmsf"].values["A"]]
    )


def test_core_values_keep_the_identity_and_follow_core_and_regions(configs, tmp_path) -> None:
    """core_rms_deviation^2 = core_rmsf^2 + core_offset^2 per replicate, over the chosen residues."""
    from polyzymd.analyses import analyze

    settings = {"selection": "all", "alignment_selection": "all", "reference_mode": "average"}
    settings |= {"core": "resid 1 2", "regions": {"tip": "resid 3"}}
    options = {"equilibration": EQUILIBRATION, "output_dir": tmp_path, "plots": False}
    got = {
        run: analyze("rmsf", [configs["A"]], run=run, settings=settings, **options)
        for run in ("core_rms_deviation", "core_rmsf", "core_offset", "tip_rmsf", "tip_offset")
    }
    values = {run: np.array(r.conditions[0].replicate_values) for run, r in got.items()}
    identity = values["core_rms_deviation"] ** 2 - values["core_rmsf"] ** 2
    assert np.max(np.abs(identity - values["core_offset"] ** 2)) < 1e-6
    assert got["tip_rmsf"].provenance.settings["residues"] == {"core": [1, 2], "tip": [3]}
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    rows = study.per_replicate(
        functions.rms_decomposition,
        pz.select("all"),
        pz.select("all"),
        pz.reference("average", "(all) or (all)", alignment="all"),
        unit="A",
        labels=lambda u: u.residues.resids,
        output_dir=tmp_path / "check",
        parts=functions.RMS_PARTS + functions.MS_PARTS,
    )
    msf = np.array(rows["msf"].values["A"])
    assert values["core_rmsf"] == pytest.approx(np.sqrt(msf[:, :2].mean(axis=1)))
    assert values["tip_rmsf"] == pytest.approx(np.sqrt(msf[:, 2]))
    with pytest.raises(ProtocolError, match="picks no residues"):
        analyze("rmsf", [configs["A"]], settings={**settings, "core": "resid 9"}, **options)
    with pytest.raises(ProtocolError, match="regions must map"):
        analyze(
            "rmsf", [configs["A"]], settings={**settings, "regions": {"core": "all"}}, **options
        )


def test_comparison_yaml_rmsf_block_is_retired(tmp_path) -> None:
    from polyzymd.config.comparison import PlotSettings, PluginSettingsContainer

    with pytest.warns(
        UserWarning, match="study.per_replicate with polyzymd.analyses.functions.rmsf"
    ):
        PluginSettingsContainer(rmsf={"selection": "name CA"})
    with pytest.warns(UserWarning, match="plot_settings.rmsf block, which is ignored"):
        PlotSettings(rmsf={"highlight_residues": [1]})


def test_bounds_warn_per_label_and_carry_to_the_mean(study, tmp_path) -> None:
    """An interval below the lower bound of one label is flagged at that label only."""
    from polyzymd.analyses.timeseries import ReplicateValues

    profile = study.per_replicate(
        _resid_value,
        pz.universe(),
        unit="A",
        labels=lambda u: u.residues.resids,
        output_dir=tmp_path,
        bounds=(0.0, None),
    )
    assert profile.bounds == (0.0, None) and profile.over_labels().bounds == (0.0, None)
    rows = {
        "A": [
            (1, np.array([0.01, 5.0]), None, None, 7, None, None),
            (2, np.array([0.5, 5.2]), None, None, 7, None, None),
        ]
    }
    near_zero = ReplicateValues(profile.source, "rmsf", "A", False, rows, [1, 2])
    near_zero.bounds = (0.0, None)
    warnings = near_zero.summary(conditions=["A"]).warnings
    assert [text for text in warnings if "extends past the bounds" in text] == [
        "the 95 percent interval of condition A at 1 extends past the bounds 0 to inf of rmsf, "
        "where a t interval is not reliable"
    ]
