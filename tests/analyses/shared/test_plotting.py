"""Tests for shared plotting helpers."""

from __future__ import annotations

from unittest.mock import patch

import matplotlib
import numpy as np
import pytest

from polyzymd.analyses.mda import (
    ArtifactStore,
    ArtifactStoreError,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.shared.plotting import (
    finite_numeric_values,
    grouped_bars,
    has_replicate_uncertainty,
    load_canonical_plot_artifacts,
    scatter_replicate_values,
    scatter_stacked_segment_replicates,
    suppress_singleton_errors,
)
from polyzymd.config.comparison import PlotSettings, PlotTheme

matplotlib.use("Agg")


def test_load_canonical_plot_artifacts_reads_configured_artifacts_only(tmp_path) -> None:
    """Plot artifact loading should ignore legacy JSON and extra run directories."""
    analysis_dir = tmp_path / "condition" / "rmsd"
    aggregated_dir = analysis_dir / "aggregated"
    run_1 = analysis_dir / "run_1"
    run_99 = analysis_dir / "run_99"
    run_99.mkdir(parents=True)

    ArtifactStore(aggregated_dir).write_condition_result(
        ConditionArtifact(
            analysis_name="rmsd",
            condition_label="condition",
            replicates=[1],
            payload={"metric": 1.0},
        )
    )
    ArtifactStore(run_1).write_replicate_result(
        ReplicateArtifact(
            analysis_name="rmsd",
            condition_label="condition",
            replicate=1,
            payload={"metric": 1.1},
        )
    )
    (run_99 / "result.json").write_text('{"artifact_type": "replicate"}', encoding="utf-8")
    (analysis_dir / "legacy_plot.json").write_text('{"metric": 99}', encoding="utf-8")

    loaded = load_canonical_plot_artifacts(analysis_dir, [1])

    assert loaded.condition_artifact is not None
    assert loaded.condition_artifact.payload["metric"] == pytest.approx(1.0)
    assert set(loaded.replicate_artifacts) == {1}
    assert loaded.replicate_artifacts[1].payload["metric"] == pytest.approx(1.1)


def test_load_canonical_plot_artifacts_rejects_legacy_result_json(tmp_path) -> None:
    """Legacy JSON at the canonical path should fail artifact validation."""
    analysis_dir = tmp_path / "condition" / "distances"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)
    (aggregated_dir / "result.json").write_text('{"pair_results": []}', encoding="utf-8")

    with pytest.raises(ArtifactStoreError, match="condition artifact"):
        load_canonical_plot_artifacts(analysis_dir, [])


def test_load_canonical_plot_artifacts_rejects_corrupt_result_json(tmp_path) -> None:
    """Corrupt canonical JSON should fail before plotters see payloads."""
    analysis_dir = tmp_path / "condition" / "sasa"
    run_dir = analysis_dir / "run_1"
    run_dir.mkdir(parents=True)
    (run_dir / "result.json").write_text("{not-json", encoding="utf-8")

    with pytest.raises(ArtifactStoreError, match="replicate artifact"):
        load_canonical_plot_artifacts(analysis_dir, [1])


def test_artifact_store_load_npz_sidecar_rejects_tampering(tmp_path) -> None:
    """NPZ sidecar loading should validate the sidecar before opening it."""
    store = ArtifactStore(tmp_path)
    sidecar = store.write_npz_sidecar(
        "sidecars/data.npz",
        values=np.asarray([1.0, 2.0], dtype=np.float64),
    )
    store.resolve_sidecar(sidecar).write_bytes(b"tampered")

    with pytest.raises(ArtifactStoreError, match="Sidecar .* mismatch"):
        store.load_npz_sidecar(sidecar)


def test_finite_numeric_values_skips_invalid_entries() -> None:
    """Finite filtering should retain numeric values and skip invalid ones."""
    values = finite_numeric_values([1.0, "bad", float("nan"), "2.5", float("inf")])

    np.testing.assert_allclose(values, [1.0, 2.5])


def test_has_replicate_uncertainty_requires_two_replicates() -> None:
    """Replicate uncertainty should require at least two finite values."""
    assert has_replicate_uncertainty([1.0]) is False
    assert has_replicate_uncertainty([1.0, 1.2]) is True
    assert has_replicate_uncertainty(n_replicates=1) is False
    assert has_replicate_uncertainty(n_replicates=2) is True


def test_suppress_singleton_errors_returns_none_for_singletons() -> None:
    """Singleton-only bars should omit uncertainty bars entirely."""
    assert suppress_singleton_errors([0.1, 0.2], [[1.0], [2.0]]) is None


def test_suppress_singleton_errors_zeroes_mixed_singletons() -> None:
    """Mixed replicate bars should keep only displayable uncertainties."""
    assert suppress_singleton_errors([0.1, 0.2], [[1.0], [2.0, 2.2]]) == [0.0, 0.2]


def test_scatter_replicate_values_vertical_filters_and_jitters() -> None:
    """Vertical scatter should jitter x positions and place replicates on y."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings(theme=PlotTheme(dot_color="red", dot_size=12, dot_alpha=0.5))
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_replicate_values(
            ax,
            [2.0],
            [[1.0, "bad", 3.0]],
            plot_settings,
            orientation="vertical",
            bar_width=0.8,
        )

    assert n_calls == 1
    np.testing.assert_allclose(mock_scatter.call_args.args[1], [1.8, 2.2])
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [1.0, 3.0])
    assert mock_scatter.call_args.kwargs["color"] == "red"
    assert mock_scatter.call_args.kwargs["s"] == 12
    assert mock_scatter.call_args.kwargs["alpha"] == 0.5
    plt.close(fig)


def test_scatter_replicate_values_horizontal_uses_values_on_x() -> None:
    """Horizontal scatter should place replicate values on x and jitter y."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        scatter_replicate_values(
            ax,
            [1.0],
            [[0.5, 1.5]],
            plot_settings,
            orientation="horizontal",
            bar_width=0.4,
        )

    np.testing.assert_allclose(mock_scatter.call_args.args[1], [0.5, 1.5])
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [0.9, 1.1])
    plt.close(fig)


def test_scatter_replicate_values_skips_disabled_theme() -> None:
    """Dot overlays should be skipped when the theme disables dot markers."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings(theme=PlotTheme(dot_size=0, dot_alpha=0.7))
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_replicate_values(ax, [0.0], [[1.0]], plot_settings)

    assert n_calls == 0
    mock_scatter.assert_not_called()
    plt.close(fig)


def test_scatter_replicate_values_rejects_length_mismatch() -> None:
    """Dot overlays should fail fast when values do not align to bars."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with pytest.raises(ValueError, match="replicate_values length must match"):
        scatter_replicate_values(ax, [0.0], [[1.0], [2.0]], plot_settings)

    plt.close(fig)


def test_scatter_stacked_segment_replicates_centers_positive_values() -> None:
    """Stacked overlays should center positive segment replicate values."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_stacked_segment_replicates(
            ax,
            1.0,
            0.4,
            [0.2, 0.4],
            plot_settings,
            bar_width=0.8,
        )

    assert n_calls == 1
    np.testing.assert_allclose(mock_scatter.call_args.args[1], [0.8, 1.2])
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [0.5, 0.6])
    plt.close(fig)


def test_scatter_stacked_segment_replicates_centers_negative_values() -> None:
    """Stacked overlays should center negative segment replicate values."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_stacked_segment_replicates(
            ax,
            1.0,
            -0.4,
            [-0.2, -0.4],
            plot_settings,
            bar_width=0.8,
        )

    assert n_calls == 1
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [-0.5, -0.6])
    plt.close(fig)


def test_scatter_stacked_segment_replicates_uses_replicate_bases() -> None:
    """Stacked overlays should use per-replicate cumulative bases."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_stacked_segment_replicates(
            ax,
            1.0,
            0.4,
            [0.2, 0.4],
            plot_settings,
            replicate_base_values=[1.0, 2.0],
            bar_width=0.8,
        )

    assert n_calls == 1
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [1.1, 2.2])
    plt.close(fig)


def test_scatter_stacked_segment_replicates_can_place_dots_at_segment_ends() -> None:
    """Stacked overlays should optionally place dots at segment ends."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_stacked_segment_replicates(
            ax,
            1.0,
            0.4,
            [0.2, 0.4],
            plot_settings,
            replicate_base_values=[1.0, 2.0],
            placement="end",
            bar_width=0.8,
        )

    assert n_calls == 1
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [1.2, 2.4])
    plt.close(fig)


def test_scatter_stacked_segment_replicates_uses_signed_replicate_bases() -> None:
    """Signed stacked overlays should choose bases from each replicate sign."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("matplotlib.axes.Axes.scatter", autospec=True) as mock_scatter:
        n_calls = scatter_stacked_segment_replicates(
            ax,
            1.0,
            0.0,
            [2.0, -4.0],
            plot_settings,
            positive_base_values=[3.0, 30.0],
            negative_base_values=[-30.0, -5.0],
            bar_width=0.8,
        )

    assert n_calls == 1
    np.testing.assert_allclose(mock_scatter.call_args.args[2], [4.0, -7.0])
    plt.close(fig)


def test_grouped_bars_uses_shared_replicate_scatter() -> None:
    """Grouped bars should delegate replicate overlays to the shared helper."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with patch("polyzymd.analyses.shared.plotting.scatter_replicate_values") as mock_scatter:
        grouped_bars(
            ax,
            np.array([0.0, 1.0]),
            [("A", [1.0, 2.0], [0.1, 0.2])],
            ["blue"],
            plot_settings,
            reference_line=None,
            replicate_values=[[[0.9, 1.1], [1.8, 2.2]]],
        )

    mock_scatter.assert_called_once()
    plt.close(fig)


def test_grouped_bars_rejects_replicate_series_mismatch() -> None:
    """Grouped bars should fail fast when replicate overlays omit a series."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with pytest.raises(ValueError, match="replicate_values length must match series length"):
        grouped_bars(
            ax,
            np.array([0.0]),
            [("A", [1.0], [0.1]), ("B", [2.0], [0.2])],
            ["blue", "orange"],
            plot_settings,
            reference_line=None,
            replicate_values=[[[0.9, 1.1]]],
        )

    plt.close(fig)


def test_grouped_bars_rejects_replicate_group_mismatch() -> None:
    """Grouped bars should validate each replicate series against x groups."""
    import matplotlib.pyplot as plt

    plot_settings = PlotSettings()
    fig, ax = plt.subplots()

    with pytest.raises(ValueError, match="replicate_values entries must match x length"):
        grouped_bars(
            ax,
            np.array([0.0, 1.0]),
            [("A", [1.0, 2.0], [0.1, 0.2])],
            ["blue"],
            plot_settings,
            reference_line=None,
            replicate_values=[[[0.9, 1.1]]],
        )

    plt.close(fig)
