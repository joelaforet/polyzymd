"""Tests for contact_events, restricted_mean_lifetime, contact_lifetimes and the lifetime runs.

The helpers are checked on hand-built masks and durations, and the
restricted mean against a product-limit Kaplan-Meier estimate written here in
numpy. contact_lifetimes and ``polyzymd analyze contacts`` are checked on the
OpenMM run directories of the systems of
tests/analyses/test_residue_contacts.py (distance: four protein residues, an
SBM and an EGM atom placed next to one residue or far away on each frame) and
tests/analyses/test_residue_occlusion.py (occlusion: an SBM and an EGM shell
that engulf one residue or sit far away). Expected events come from the
schedule of which residue each polymer residue touches on each frame,
through :func:`_expected_events` below, which never calls the code under
test.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import (
    LIFETIME_PARTS,
    contact_events,
    contact_lifetimes,
    restricted_mean_lifetime,
)
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config
from tests.analyses import test_residue_contacts as distance_system
from tests.analyses import test_residue_occlusion as occlusion_system

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0ns"
MEAN, EVENTS, CENSORED = range(3)
TYPES = ("EGM", "SBM")


# ---------------------------------------------------------------------------
# An independent product-limit estimate and event finder
# ---------------------------------------------------------------------------


def _product_limit_area(durations, censored, horizon: float) -> float:
    """Area to ``horizon`` under S(t) = prod over event times t_i <= t of (1 - d_i / n_i).

    ``d_i`` is the number of uncensored durations equal to ``t_i`` and ``n_i``
    the number of durations, censored or not, of at least ``t_i``.
    """
    durations = np.asarray(durations, dtype=float)
    censored = np.asarray(censored, dtype=bool)
    times = np.unique(durations[~censored])
    level, area, previous = 1.0, 0.0, 0.0
    for t in times:
        if t >= horizon:
            break
        area += level * (t - previous)
        n_at_risk = np.sum(durations >= t)
        n_ended = np.sum(durations[~censored] == t)
        level *= 1.0 - n_ended / n_at_risk
        previous = t
    return float(area + level * (horizon - previous))


def _runs(present: list[bool], gap: int) -> list[tuple[int, bool]]:
    """Runs of True in ``present`` as (length, censored), after filling interior absences <= ``gap``."""
    present = list(present)
    n = len(present)
    on = [i for i, value in enumerate(present) if value]
    for a, b in zip(on, on[1:]):
        if 1 < b - a <= gap + 1:
            present[a:b] = [True] * (b - a)
    runs, start = [], None
    for i, value in enumerate([*present, False]):
        if value and start is None:
            start = i
        elif not value and start is not None:
            runs.append((i - start, start == 0 or i == n))
            start = None
    return runs


def _expected_events(schedule, residues, gap: int = 0) -> dict[str, list[tuple[int, bool]]]:
    """Events of the polymer, EGM and SBM on each residue of ``residues``, pooled."""
    events = {"polymer": [], "EGM": [], "SBM": []}
    for r in residues:
        events["polymer"] += _runs([r in (sbm, egm) for sbm, egm in schedule], gap)
        events["SBM"] += _runs([sbm == r for sbm, _ in schedule], gap)
        events["EGM"] += _runs([egm == r for _, egm in schedule], gap)
    return events


def _expected_row(events, step_ns: float, n_frames: int) -> list[float]:
    """mean_lifetime, n_events and censored_fraction of ``events`` at ``step_ns`` per frame."""
    if not events:
        return [float("nan"), 0.0, float("nan")]
    lengths = np.array([length for length, _ in events], dtype=float)
    censored = np.array([flag for _, flag in events])
    area = _product_limit_area(lengths * step_ns, censored, n_frames * step_ns)
    return [area, float(len(events)), float(np.mean(censored))]


def _expected_table(schedule, residues, step_ns, gap=0, frames=None) -> np.ndarray:
    """The (3, 3) table of contact_lifetimes with types EGM and SBM, from the schedule."""
    events = _expected_events(schedule, residues, gap)
    return np.array(
        [_expected_row(events[name], step_ns, len(schedule)) for name in ("polymer", *TYPES)]
    ).T


# ---------------------------------------------------------------------------
# contact_events
# ---------------------------------------------------------------------------


def _mask(*columns: str) -> np.ndarray:
    """A frames x columns mask from strings such as '.XX.', one per column."""
    return np.array([[c == "X" for c in column] for column in columns]).T


def test_contact_events_finds_runs_and_their_lengths() -> None:
    lengths, censored = contact_events(_mask(".XX..XXX."))

    assert lengths.tolist() == [2, 3]
    assert censored.tolist() == [False, False]


@pytest.mark.parametrize(
    ("column", "lengths", "censored"),
    [
        ("XX...", [2], [True]),
        ("...XX", [2], [True]),
        ("XXXXX", [5], [True]),
        ("X.X.X", [1, 1, 1], [True, False, True]),
        (".X.X.", [1, 1], [False, False]),
    ],
)
def test_contact_events_censors_runs_at_the_first_or_last_frame(column, lengths, censored) -> None:
    found, flags = contact_events(_mask(column))

    assert found.tolist() == lengths
    assert flags.tolist() == censored


def test_contact_events_pools_the_columns_in_order() -> None:
    lengths, censored = contact_events(_mask("XX.X.", "....."), gap=0)
    assert (lengths.tolist(), censored.tolist()) == ([2, 1], [True, False])

    lengths, censored = contact_events(_mask(".XX..", "X..XX", "....."))

    assert lengths.tolist() == [2, 1, 2]
    assert censored.tolist() == [False, True, True]
    assert lengths.dtype.kind == "i" and censored.dtype == bool


def test_contact_events_of_an_all_false_mask_are_empty() -> None:
    lengths, censored = contact_events(_mask(".....", "....."), gap=2)

    assert lengths.tolist() == [] and censored.tolist() == []
    assert lengths.dtype.kind == "i" and censored.dtype == bool


@pytest.mark.parametrize(
    ("column", "gap", "lengths", "censored"),
    [
        (".X.X.", 1, [3], [False]),
        (".X..X.", 1, [1, 1], [False, False]),
        (".X..X.", 2, [4], [False]),
        (".X...X.", 2, [1, 1], [False, False]),
        (".X...X.", 3, [5], [False]),
        (".X.X..X.", 1, [3, 1], [False, False]),
        (".X.X..X.", 2, [6], [False]),
        ("X.X.", 1, [3], [True]),
    ],
)
def test_contact_events_fill_interior_absences_up_to_gap(column, gap, lengths, censored) -> None:
    found, flags = contact_events(_mask(column), gap=gap)

    assert found.tolist() == lengths
    assert flags.tolist() == censored


@pytest.mark.parametrize("column", [".XX.", "..XX", "XX..", "..X.."])
def test_contact_events_do_not_fill_absences_at_the_ends(column) -> None:
    unfilled = contact_events(_mask(column))
    filled = contact_events(_mask(column), gap=3)

    assert filled[0].tolist() == unfilled[0].tolist()
    assert filled[1].tolist() == unfilled[1].tolist()


def test_contact_events_fill_each_column_on_its_own() -> None:
    """A presence in one column does not bridge a gap in another."""
    lengths, _ = contact_events(_mask("X.X..", ".X..X"), gap=1)

    assert lengths.tolist() == [3, 1, 1]


def test_contact_events_match_the_independent_finder_on_random_masks() -> None:
    rng = np.random.default_rng(7)
    for gap in (0, 1, 2, 4):
        mask = rng.random((40, 6)) < 0.45
        lengths, censored = contact_events(mask, gap=gap)
        expected = [run for column in mask.T for run in _runs(column.tolist(), gap)]
        assert lengths.tolist() == [length for length, _ in expected], gap
        assert censored.tolist() == [flag for _, flag in expected], gap


# ---------------------------------------------------------------------------
# restricted_mean_lifetime
# ---------------------------------------------------------------------------


def test_without_censoring_and_a_long_horizon_it_is_the_mean() -> None:
    durations = [0.3, 1.2, 0.5, 0.5, 2.0]

    for horizon in (2.0, 5.0, 100.0):
        value = restricted_mean_lifetime(durations, [False] * 5, horizon)
        assert value == pytest.approx(np.mean(durations), rel=1e-12)


@pytest.mark.parametrize("horizon", [0.1, 0.5, 0.8, 1.5])
def test_without_censoring_it_is_the_mean_of_durations_cut_at_the_horizon(horizon) -> None:
    durations = np.array([0.3, 1.2, 0.5, 0.5, 2.0])

    value = restricted_mean_lifetime(durations, np.zeros(5, dtype=bool), horizon)

    assert value == pytest.approx(np.mean(np.minimum(durations, horizon)), rel=1e-12)


@pytest.mark.parametrize("seed", range(6))
def test_it_matches_an_independent_product_limit_estimate(seed) -> None:
    """Random integer durations give ties, among events and between events and censored ones."""
    rng = np.random.default_rng(seed)
    n = int(rng.integers(5, 60))
    durations = rng.integers(1, 12, size=n) * 0.04
    censored = rng.random(n) < 0.35
    censored[0] = False  # at least one event

    for horizon in (0.1, 0.25, float(durations.max()), 0.6):
        value = restricted_mean_lifetime(durations, censored, horizon)
        assert value == pytest.approx(_product_limit_area(durations, censored, horizon), rel=1e-10)


def test_a_censored_duration_tied_with_an_event_is_still_at_risk() -> None:
    """At 0.2 three durations are at risk, one ends: S = 3/4 * 2/3 after 0.2."""
    value = restricted_mean_lifetime([0.2, 0.1, 0.2, 0.2], [True, False, False, True], 0.8)

    assert value == pytest.approx(0.1 + 0.1 * 3 / 4 + 0.6 * 1 / 2, rel=1e-12)


def test_it_is_nan_without_durations() -> None:
    assert np.isnan(restricted_mean_lifetime([], [], 1.0))
    assert np.isnan(restricted_mean_lifetime(np.array([]), np.array([], dtype=bool), 1.0))


def test_a_censored_longest_duration_keeps_the_tail_flat_to_the_horizon() -> None:
    """After 0.1 the survival stays at 1/2, since the one longer duration is censored."""
    value = restricted_mean_lifetime([0.1, 0.3], [False, True], 1.0)

    assert value == pytest.approx(0.1 + 0.9 * 0.5, rel=1e-12)
    assert restricted_mean_lifetime([0.1, 0.3], [False, False], 1.0) == pytest.approx(0.2)


def test_only_censored_durations_give_the_horizon() -> None:
    assert restricted_mean_lifetime([0.2, 0.1], [True, True], 0.4) == pytest.approx(0.4)


# ---------------------------------------------------------------------------
# contact_lifetimes on OpenMM run directories
# ---------------------------------------------------------------------------

#: Eight frames of (SBM residue, EGM residue), residues 1 to 4 of the distance
#: system. Polymer events: residue 1 on frames 0-1 (censored) and 3, residue 2
#: on 4-5 and residue 3 on 6-7 (censored); residue 4 is never touched. SBM:
#: residue 1 on 0-1 (censored) and 3, residue 2 on 4-5. EGM: residue 2 on 4
#: and residue 3 on 6-7 (censored).
SCHEDULE = [(1, None), (1, None), (None, None), (1, None), (2, 2), (2, None), (None, 3), (None, 3)]
RESIDUES = [1, 2, 3, 4]


def _write_distance(config: Path, replicate: int, schedule, **kwargs) -> None:
    write_openmm_frames(
        config,
        replicate,
        distance_system._frames(schedule),
        distance_system.RESINDEX,
        names=distance_system.NAMES,
        resnames=distance_system.RESNAMES,
        elements=distance_system.ELEMENTS,
        chain_ids=distance_system.CHAIN_IDS,
        **kwargs,
    )


def _universe(root: Path, schedule, *, dt_ps: float = 100.0, occlusion: bool = False):
    """The universe of one run directory of ``schedule``, loaded as a Study loads it."""
    config = write_simulation_config(root, scratch=root / "scratch")
    if occlusion:
        occlusion_system._write(config, 1, schedule, bonds=occlusion_system.BONDS, dt_ps=dt_ps)
    else:
        _write_distance(config, 1, schedule, dt_ps=dt_ps)
    study = pz.Study.from_configs({"x": config}, equilibration=EQUILIBRATION)
    return next(iter(study)).replicates[0].universe()


def _groups(universe):
    return universe.select_atoms("chainid A"), universe.select_atoms("chainid C")


def _lifetimes(universe, frames=None, **kwargs) -> np.ndarray:
    protein, polymer = _groups(universe)
    frames = range(universe.trajectory.n_frames) if frames is None else frames
    return contact_lifetimes(protein, polymer, frames, **kwargs)


@pytest.fixture(scope="module")
def distance_universe(tmp_path_factory):
    return _universe(tmp_path_factory.mktemp("distance"), SCHEDULE)


def test_distance_lifetimes_equal_the_hand_computed_values(distance_universe) -> None:
    result = _lifetimes(distance_universe, method="distance", types=TYPES)

    assert result.shape == (len(LIFETIME_PARTS), 3)
    assert LIFETIME_PARTS == ("mean_lifetime", "n_events", "censored_fraction")
    # Polymer: 0.2 (c), 0.1, 0.2, 0.2 (c) ns up to 0.8 ns; S = 1, then 3/4 after
    # 0.1, then 1/2 after 0.2.
    assert result[:, 0] == pytest.approx([0.1 + 0.1 * 3 / 4 + 0.6 / 2, 4, 0.5], rel=1e-6)
    # EGM: 0.1, 0.2 (c): S = 1/2 after 0.1.
    assert result[:, 1] == pytest.approx([0.1 + 0.7 / 2, 2, 0.5], rel=1e-6)
    # SBM: 0.2 (c), 0.1, 0.2: S = 2/3 after 0.1, 1/3 after 0.2.
    assert result[:, 2] == pytest.approx([0.1 + 0.1 * 2 / 3 + 0.6 / 3, 3, 1 / 3], rel=1e-6)
    expected = _expected_table(SCHEDULE, RESIDUES, 0.1)
    assert result == pytest.approx(expected, rel=1e-6)


def test_without_types_there_is_one_column_and_an_absent_type_has_no_event(
    distance_universe,
) -> None:
    plain = _lifetimes(distance_universe, method="distance")
    absent = _lifetimes(distance_universe, method="distance", types=("XYZ",))

    assert plain.shape == (3, 1)
    assert absent[:, 0] == pytest.approx(plain[:, 0])
    assert np.isnan(absent[MEAN, 1]) and np.isnan(absent[CENSORED, 1])
    assert absent[EVENTS, 1] == 0


def test_tolerance_of_one_frame_spacing_fills_the_one_frame_gap(distance_universe) -> None:
    """Residue 1's frames 0-1 and 3 become one censored event of 4 frames."""
    filled = _lifetimes(distance_universe, method="distance", types=TYPES, tolerance_ps=100.0)
    short = _lifetimes(distance_universe, method="distance", types=TYPES, tolerance_ps=99.0)
    plain = _lifetimes(distance_universe, method="distance", types=TYPES)

    # Polymer: 0.4 (c), 0.2, 0.2 (c) ns: S = 2/3 after 0.2.
    assert filled[:, 0] == pytest.approx([0.2 + 0.6 * 2 / 3, 3, 2 / 3], rel=1e-6)
    assert filled == pytest.approx(_expected_table(SCHEDULE, RESIDUES, 0.1, gap=1), rel=1e-6)
    assert short == pytest.approx(plain, rel=1e-12)
    assert _lifetimes(distance_universe, method="distance", types=TYPES, tolerance_ps=0.0) == (
        pytest.approx(plain, rel=1e-12)
    )


def test_a_tolerance_equal_to_a_spacing_stored_just_above_it_fills_one_frame(tmp_path) -> None:
    """A DCD written at 30 ps stores a spacing a little above 30 ps, and 30 ps still fills one frame."""
    universe = _universe(tmp_path, SCHEDULE, dt_ps=30.0)
    times = np.array([ts.time for ts in universe.trajectory])
    assert np.median(np.diff(times)) > 30.0  # the case the rounding allowance is for

    filled = _lifetimes(universe, method="distance", types=TYPES, tolerance_ps=30.0)
    unfilled = _lifetimes(universe, method="distance", types=TYPES, tolerance_ps=29.0)

    assert filled == pytest.approx(_expected_table(SCHEDULE, RESIDUES, 0.03, gap=1), rel=1e-6)
    assert unfilled == pytest.approx(_expected_table(SCHEDULE, RESIDUES, 0.03), rel=1e-6)
    assert (filled[EVENTS, 0], unfilled[EVENTS, 0]) == (3, 4)


def test_every_other_frame_doubles_the_spacing(distance_universe) -> None:
    result = _lifetimes(distance_universe, frames=[0, 2, 4, 6], method="distance", types=TYPES)

    # Polymer: residue 1 on frame 0 (c), 2 on frame 4 and 3 on frame 6 (c), 0.2 ns each.
    assert result[:, 0] == pytest.approx([0.2 + 0.6 * 2 / 3, 3, 2 / 3], rel=1e-6)
    assert result == pytest.approx(_expected_table(SCHEDULE[::2], RESIDUES, 0.2), rel=1e-6)


def test_distance_options_reach_the_contacts(distance_universe) -> None:
    """The atoms sit 3.04 Å apart, so a 3 Å cutoff finds no contact."""
    result = _lifetimes(distance_universe, method="distance", types=TYPES, cutoff=3.0, pbc=False)

    assert result[EVENTS].tolist() == [0, 0, 0]
    assert np.isnan(result[MEAN]).all()


#: Four frames of the occlusion system, measured residues 2 to 5. Polymer:
#: residue 2 on 0-1 (censored), 3 on 1-2, 4 on 3 (censored). SBM: residue 2
#: on 0-1 and 4 on 3, both censored. EGM: residue 3 on 1-2.
OCCLUSION_SCHEDULE = [(2, None), (2, 3), (None, 3), (4, None)]


def test_occlusion_lifetimes_equal_the_hand_computed_values(tmp_path) -> None:
    universe = _universe(tmp_path, OCCLUSION_SCHEDULE, occlusion=True)

    result = _lifetimes(universe, types=TYPES)
    distance = _lifetimes(universe, method="distance", types=TYPES, cutoff=1.0)

    # Polymer: 0.2 (c), 0.2, 0.1 (c) ns up to 0.4 ns: S = 1/2 after 0.2.
    assert result[:, 0] == pytest.approx([0.2 + 0.2 / 2, 3, 2 / 3], rel=1e-6)
    # EGM: one uncensored 0.2 ns event.
    assert result[:, 1] == pytest.approx([0.2, 1, 0.0], rel=1e-6)
    # SBM: both events censored, so S = 1 to the horizon.
    assert result[:, 2] == pytest.approx([0.4, 2, 1.0], rel=1e-6)
    expected = _expected_table(OCCLUSION_SCHEDULE, occlusion_system.MEASURED, 0.1)
    assert result == pytest.approx(expected, rel=1e-6)
    assert distance[EVENTS].tolist() == [0, 0, 0]  # the shells are 2.6 Å away or more


def test_occlusion_options_reach_the_contacts(tmp_path) -> None:
    """A threshold of 0 needs a residue to lose more than all its area, which never happens."""
    universe = _universe(tmp_path, OCCLUSION_SCHEDULE, occlusion=True)

    result = _lifetimes(universe, method="occlusion", types=TYPES, threshold=0.0)

    assert result[EVENTS].tolist() == [0, 0, 0]


def test_uneven_frame_spacing_is_refused(distance_universe) -> None:
    with pytest.raises(ProtocolError, match="not evenly spaced") as info:
        _lifetimes(distance_universe, frames=[0, 1, 3], method="distance")
    assert info.value.hint


@pytest.mark.parametrize("frames", [[], [2]])
def test_fewer_than_two_frames_are_refused(distance_universe, frames) -> None:
    with pytest.raises(ProtocolError, match="at least two") as info:
        _lifetimes(distance_universe, frames=frames, method="distance")
    assert info.value.hint


def test_an_unknown_method_is_refused(distance_universe) -> None:
    with pytest.raises(ProtocolError, match="method must be 'occlusion' or 'distance'"):
        _lifetimes(distance_universe, method="sasa")


@pytest.mark.parametrize(
    ("method", "option", "match"),
    [
        ("distance", {"threshold": 0.3}, "method distance takes no option threshold"),
        ("distance", {"max_asa": "empirical"}, "method distance takes no option max_asa"),
        ("occlusion", {"cutoff": 4.0}, "method occlusion takes no option cutoff"),
        ("occlusion", {"bogus": 1}, "method occlusion takes no option bogus"),
    ],
)
def test_an_option_of_the_other_method_is_refused(distance_universe, method, option, match):
    with pytest.raises(ProtocolError, match=match) as info:
        _lifetimes(distance_universe, method=method, **option)
    assert info.value.hint


# ---------------------------------------------------------------------------
# polyzymd analyze contacts: the lifetime runs
# ---------------------------------------------------------------------------

#: Three replicates of eight frames. Replicate 1 is SCHEDULE; in 2 SBM holds
#: residue 1 throughout; in 3 SBM hops between residues every other frame.
#: EGM touches a residue only in replicate 1.
REPLICATES = {
    1: SCHEDULE,
    2: [(1, None)] * 8,
    3: [(1, None), (None, None), (2, None), (None, None)] * 2,
}
LIFETIME_RUNS = {
    "mean_lifetime": (MEAN, "polymer"),
    "EGM_mean_lifetime": (MEAN, "EGM"),
    "SBM_mean_lifetime": (MEAN, "SBM"),
    "lifetime_events": (EVENTS, "polymer"),
    "censored_fraction": (CENSORED, "polymer"),
}


@pytest.fixture(scope="module")
def config(tmp_path_factory) -> Path:
    root = tmp_path_factory.mktemp("lifetimes")
    path = write_simulation_config(root, scratch=root / "scratch")
    for replicate, schedule in REPLICATES.items():
        _write_distance(path, replicate, schedule)
    return path


def _options(tmp_path: Path, settings: dict | None = None, **extra):
    return {
        "equilibration": EQUILIBRATION,
        "output_dir": tmp_path,
        "plots": False,
        "labels": ["A"],
        "settings": {"method": "distance", **(settings or {})},
        **extra,
    }


def _expected_values(run: str, gap: int = 0) -> list[float]:
    row, group = LIFETIME_RUNS[run]
    column = ("polymer", *TYPES).index(group)
    return [
        float(_expected_table(REPLICATES[r], RESIDUES, 0.1, gap=gap)[row, column])
        for r in sorted(REPLICATES)
    ]


@pytest.mark.parametrize("run", list(LIFETIME_RUNS))
def test_analyze_lifetime_runs_equal_the_hand_computed_values(config, tmp_path, run) -> None:
    report = analyze("contacts", [config], run=run, **_options(tmp_path))

    assert (report.analysis, report.run) == ("contacts", run)
    assert report.all_runs[-5:] == list(LIFETIME_RUNS)
    assert report.unit == ("ns" if run.endswith("mean_lifetime") else None)
    (row,) = report.conditions
    expected = _expected_values(run)
    assert row.replicates == [1, 2, 3]
    assert np.allclose(row.replicate_values, expected, rtol=1e-6, equal_nan=True)


def test_analyze_lifetime_values_by_hand_for_the_three_replicates(config, tmp_path) -> None:
    """Replicate 2: one censored 0.8 ns event; replicate 3: four 0.1 ns events, one censored."""
    lifetime = analyze("contacts", [config], run="mean_lifetime", **_options(tmp_path))
    events = analyze("contacts", [config], run="lifetime_events", **_options(tmp_path))
    censored = analyze("contacts", [config], run="censored_fraction", **_options(tmp_path))

    # Replicate 3: three of four events end at 0.1 ns, so S = 1/4 after 0.1 ns.
    assert lifetime.conditions[0].replicate_values == pytest.approx(
        [0.475, 0.8, 0.1 + 0.7 / 4], rel=1e-6
    )
    assert events.conditions[0].replicate_values == [4.0, 1.0, 4.0]
    assert censored.conditions[0].replicate_values == pytest.approx([0.5, 1.0, 0.25])


def test_analyze_bounds_flag_intervals_past_zero_and_one(config, tmp_path) -> None:
    """Censored fractions 0.5, 1 and 0.25 give a t interval past 1; lifetimes one past 0."""
    censored = analyze("contacts", [config], run="censored_fraction", **_options(tmp_path))
    lifetime = analyze("contacts", [config], run="SBM_mean_lifetime", **_options(tmp_path))

    assert any("bounds 0 to 1 of censored_fraction" in text for text in censored.warnings)
    low, _ = lifetime.conditions[0].ci95
    assert low < 0
    assert any("bounds 0 to inf of SBM_mean_lifetime" in text for text in lifetime.warnings)


def test_analyze_warns_when_a_replicate_has_no_contact_event(config, tmp_path) -> None:
    report = analyze("contacts", [config], run="EGM_mean_lifetime", **_options(tmp_path))

    values = report.conditions[0].replicate_values
    assert np.isfinite(values[0]) and np.isnan(values[1]) and np.isnan(values[2])
    (warning,) = [text for text in report.warnings if "no contact event" in text]
    assert "A replicate 2, A replicate 3" in warning
    assert "for EGM" in warning and "EGM_mean_lifetime is undefined (nan)" in warning
    plain = analyze("contacts", [config], run="mean_lifetime", **_options(tmp_path))
    assert not any("no contact event" in text for text in plain.warnings)


def test_analyze_tolerance_fills_gaps_and_is_in_the_provenance(config, tmp_path) -> None:
    default = analyze("contacts", [config], run="mean_lifetime", **_options(tmp_path / "0"))
    report = analyze(
        "contacts",
        [config],
        run="lifetime_events",
        **_options(tmp_path / "1", {"tolerance_ps": 100.0}),
    )

    assert default.provenance.settings["tolerance_ps"] == 0.0
    assert report.provenance.settings["tolerance_ps"] == 100.0
    assert report.provenance.settings["method"] == "distance"
    assert report.provenance.settings["polymer_types_found"] == list(TYPES)
    assert report.conditions[0].replicate_values == _expected_values("lifetime_events", gap=1)
    assert report.conditions[0].replicate_values == [3.0, 1.0, 4.0]


def test_analyze_stride_doubles_the_spacing_of_the_lifetimes(config, tmp_path) -> None:
    report = analyze("contacts", [config], run="mean_lifetime", **_options(tmp_path, stride=2))

    expected = [
        float(_expected_table(REPLICATES[r][::2], RESIDUES, 0.2)[MEAN, 0]) for r in (1, 2, 3)
    ]
    assert report.conditions[0].replicate_values == pytest.approx(expected, rel=1e-6)


def test_analyze_refuses_a_negative_tolerance(config) -> None:
    with pytest.raises(ProtocolError, match="tolerance_ps must be at least 0") as info:
        analyze(
            "contacts",
            [config],
            run="mean_lifetime",
            equilibration=EQUILIBRATION,
            settings={"method": "distance", "tolerance_ps": -1.0},
        )
    assert "tolerance_ps=0" in info.value.hint


def test_analyze_refuses_a_negative_tolerance_for_any_run(config, tmp_path) -> None:
    with pytest.raises(ProtocolError, match="tolerance_ps must be at least 0"):
        analyze("contacts", [config], **_options(tmp_path, {"tolerance_ps": -1.0}))


@pytest.mark.parametrize(
    ("method", "extra"),
    [("distance", {"method": "distance"}), ("occlusion", {})],
)
def test_study_per_replicate_at_the_defaults_reuses_what_analyze_stored(
    tmp_path, method, extra
) -> None:
    """analyze passes only non-default options, so a plain Python call finds its records."""
    config = write_simulation_config(tmp_path / "c", scratch=tmp_path / "c" / "scratch")
    for replicate in (1, 2, 3):
        if method == "distance":
            _write_distance(config, replicate, REPLICATES[replicate])
        else:
            occlusion_system._write(config, replicate, OCCLUSION_SCHEDULE)
    report = analyze(
        "contacts",
        [config],
        run="mean_lifetime",
        **{**_options(tmp_path), "settings": {"method": method}},
    )
    stored = sorted((tmp_path / "polyzymd_results" / "contact_lifetimes").rglob("*.npz"))
    assert stored
    before = [path.stat().st_mtime_ns for path in stored]
    protein = report.provenance.settings["protein_selection"]
    polymer = report.provenance.settings["polymer_selection"]
    types = report.provenance.settings["polymer_types_found"]
    study = pz.Study.from_configs({"A": config}, equilibration=EQUILIBRATION)

    rows = study.per_replicate(
        functions.contact_lifetimes,
        pz.select(protein),
        pz.select(polymer),
        unit=None,
        labels=["polymer", *types],
        name="contact_lifetimes",
        output_dir=tmp_path,
        parts=list(LIFETIME_PARTS),
        types=types,
        **extra,
    )

    assert [path.stat().st_mtime_ns for path in stored] == before
    polymer_lifetimes = [float(values[0]) for values in rows["mean_lifetime"].values["A"]]
    assert polymer_lifetimes == pytest.approx(report.conditions[0].replicate_values)


def test_cli_lifetime_draws_the_comparison_figure(config, tmp_path) -> None:
    arguments = ["contacts", "-c", str(config), "--eq", EQUILIBRATION]
    arguments += ["--output-dir", str(tmp_path), "--set", "method=distance"]

    result = CliRunner().invoke(analyze_command, [*arguments, "--run", "mean_lifetime"])

    assert result.exit_code == 0, result.output
    assert result.stdout.startswith("# polyzymd analyze contacts")
    assert {path.name for path in (tmp_path / "figures" / "contacts").iterdir()} == {
        "contacts_mean_lifetime_comparison.png"
    }
