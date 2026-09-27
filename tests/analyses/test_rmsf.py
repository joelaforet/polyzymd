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
        about_reference=mode == "external",
        unit="A",
        labels=lambda u: u.residues.resids,
        output_dir=tmp_path,
    )


def _by_hand(coordinates: np.ndarray, reference: np.ndarray, about_reference: bool) -> np.ndarray:
    """Align with AlignTraj onto ``reference`` and take RMSF, or deviations, per residue."""
    from MDAnalysis.analysis import align, rms

    mobile = mda.Universe.empty(6, n_residues=3, atom_resindex=RESINDEX, trajectory=True)
    mobile.load_new(coordinates.astype(np.float32), format="MEMORY")
    ref = mda.Merge(mobile.atoms)
    ref.load_new(reference.astype(np.float32)[np.newaxis], format="MEMORY")
    align.AlignTraj(mobile, ref, select="all", in_memory=True).run()
    if about_reference:
        moved = np.array([mobile.atoms.positions for _ in mobile.trajectory], float)
        per_atom = np.sqrt(np.mean(np.sum((moved - reference) ** 2, axis=2), axis=0))
    else:
        per_atom = rms.RMSF(mobile.atoms).run().results.rmsf
    return per_atom.reshape(3, 2).mean(axis=1)


@pytest.mark.parametrize("about_reference", [False, True])
def test_rmsf_equals_aligntraj_and_rms_rmsf(about_reference: bool) -> None:
    """rmsf gives what AlignTraj and rms.RMSF give, without moving the trajectory."""
    coordinates = _frames(3, 0.5)
    universe = mda.Universe.empty(6, n_residues=3, atom_resindex=RESINDEX, trajectory=True)
    universe.load_new(coordinates.astype(np.float32), format="MEMORY")
    before = universe.trajectory.coordinate_array.copy()
    ref = mda.Merge(universe.atoms)
    ref.load_new(coordinates[4][np.newaxis].astype(np.float32), format="MEMORY")
    frames = np.arange(2, 10)
    measured = rmsf(universe.atoms, universe.atoms[:4], ref.atoms, frames, about_reference)
    fit_only = _by_hand(coordinates[2:], coordinates[4], about_reference)
    if not about_reference:
        # Fitting on residues 1 and 2 only gives other values than fitting on all atoms.
        everything = rmsf(universe.atoms, universe.atoms, ref.atoms, frames)
        assert everything == pytest.approx(fit_only, abs=1e-5)
        assert not np.allclose(measured, everything, atol=1e-3)
    else:
        assert rmsf(universe.atoms, universe.atoms, ref.atoms, frames, True) == pytest.approx(
            fit_only, abs=1e-5
        )
    assert np.array_equal(universe.trajectory.coordinate_array, before)
    with pytest.raises(ProtocolError, match="reference has 4 atoms"):
        rmsf(universe.atoms, universe.atoms, ref.atoms[:4], frames)


def test_per_replicate_profile_is_labelled_by_residue(study, tmp_path) -> None:
    """Each replicate's profile equals the hand-aligned RMSF, and summary rows are per residue."""
    profile = _profile(study, tmp_path, "frame", frame=1)
    assert profile.labels == [1, 2, 3]
    replicate = study["B"].replicates[1]
    u = replicate.universe()
    coordinates = np.array([u.atoms.positions for _ in u.trajectory[replicate.frames]], float)
    expected = _by_hand(coordinates, coordinates[0], False)
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
    assert lines[0].startswith("# polyzymd analyze rmsf  metric mean_rmsf  unit A  run mean_rmsf")
    assert lines[-1].startswith("verdict: B larger mean_rmsf than A")
    folder = tmp_path / "figures" / "rmsf"
    assert {path.name for path in folder.iterdir()} == {"rmsf_profile.png", "rmsf_comparison.png"}
    per_residue = CliRunner().invoke(analyze_command, [*arguments, "--run", "per_residue"])
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


def test_comparison_yaml_rmsf_block_is_retired(tmp_path) -> None:
    from polyzymd.config.comparison import PlotSettings, PluginSettingsContainer

    with pytest.warns(
        UserWarning, match="study.per_replicate with polyzymd.analyses.functions.rmsf"
    ):
        PluginSettingsContainer(rmsf={"selection": "name CA"})
    with pytest.warns(UserWarning, match="plot_settings.rmsf block, which is ignored"):
        PlotSettings(rmsf={"highlight_residues": [1]})
