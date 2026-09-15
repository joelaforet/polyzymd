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
    _finite_numeric_values,
    get_condition_color_map,
    get_condition_colors,
    get_palette_colors,
    grouped_bars,
    load_canonical_plot_artifacts,
    order_condition_labels,
    scatter_replicate_values,
)
from polyzymd.config.comparison import PlotSettings, PlotTheme

matplotlib.use("Agg")


def _rgba(color):
    """Normalize matplotlib color specs for assertions.

    Parameters
    ----------
    color : Any
        Matplotlib-compatible color specification.

    Returns
    -------
    tuple[float, float, float, float]
        RGBA tuple for stable comparisons.
    """
    from matplotlib.colors import to_rgba

    return to_rgba(color)


def test_semantic_order_disabled_preserves_input_order() -> None:
    """Disabled semantic colors should preserve existing condition order."""
    plot_settings = PlotSettings(
        semantic_colors={"enabled": False, "order": ["B", "A"]},
    )

    assert order_condition_labels(["A", "B", "C"], plot_settings) == ["A", "B", "C"]


def test_semantic_colors_disabled_use_noncanonical_palette() -> None:
    """Disabled semantic colors should return palette colors."""
    labels = ["A", "B", "C"]
    plot_settings = PlotSettings(color_palette="tab10")

    assert get_condition_colors(labels, plot_settings) == get_palette_colors(3, plot_settings)


def test_invalid_global_palette_falls_back_to_tab10(caplog: pytest.LogCaptureFixture) -> None:
    """Invalid global palettes should warn and use a safe tab10 fallback."""
    plot_settings = PlotSettings(color_palette="not-a-real-palette")

    colors = get_palette_colors(3, plot_settings)

    np.testing.assert_allclose(
        [_rgba(color) for color in colors],
        [_rgba(matplotlib.colormaps["tab10"](fraction)) for fraction in (0.0, 0.5, 1.0)],
    )
    assert "not-a-real-palette" in caplog.text
    assert "Falling back to 'tab10'" in caplog.text


def test_semantic_explicit_and_condition_order() -> None:
    """Explicit order should lead, then condition order, then original order."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "order": ["C", "Missing"],
            "conditions": {
                "D": {"order": 1},
                "B": {"order": 0},
            },
        }
    )

    assert order_condition_labels(["A", "B", "C", "D", "E"], plot_settings) == [
        "C",
        "B",
        "D",
        "A",
        "E",
    ]


def test_semantic_manual_condition_and_control_precedence() -> None:
    """Manual colors should override condition and control colors."""
    labels = ["Control", "Treated", "Manual"]
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "control_color": "black",
            "manual_colors": {"Manual": "purple", "Control": "orange"},
            "conditions": {
                "Control": {"role": "control", "color": "red"},
                "Treated": {"color": "blue"},
                "Manual": {"color": "green"},
            },
        }
    )

    color_map = get_condition_color_map(labels, plot_settings, control_label="Control")

    assert color_map["Control"] == "orange"
    assert color_map["Treated"] == "blue"
    assert color_map["Manual"] == "purple"


def test_semantic_control_role_uses_control_color() -> None:
    """Control role should resolve to the configured control color."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "control_color": "black",
            "conditions": {"Reference": {"role": "Control"}},
        }
    )

    assert get_condition_color_map(["Reference"], plot_settings)["Reference"] == "black"


def test_semantic_ordinal_family_colors_use_value_order() -> None:
    """Ordinal family colors should sample configured value order."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "conditions": {
                "Low": {"family": "dose", "value": "low"},
                "High": {"family": "dose", "value": "high"},
            },
            "families": {
                "dose": {
                    "scale": "ordinal",
                    "colormap": "viridis",
                    "value_order": ["low", "medium", "high"],
                    "colormap_range": [0.0, 1.0],
                }
            },
        }
    )

    color_map = get_condition_color_map(["Low", "High"], plot_settings)

    assert _rgba(color_map["Low"]) == pytest.approx(_rgba(matplotlib.colormaps["viridis"](0.0)))
    assert _rgba(color_map["High"]) == pytest.approx(_rgba(matplotlib.colormaps["viridis"](1.0)))


def test_semantic_linear_family_colors_use_numeric_values() -> None:
    """Linear family colors should normalize numeric condition values."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "conditions": {
                "Zero": {"family": "fraction", "value": 0.0},
                "Half": {"family": "fraction", "value": 0.5},
                "One": {"family": "fraction", "value": 1.0},
            },
            "families": {
                "fraction": {
                    "scale": "linear",
                    "colormap": "plasma",
                    "vmin": 0.0,
                    "vmax": 1.0,
                    "colormap_range": [0.2, 0.8],
                }
            },
        }
    )

    color_map = get_condition_color_map(["Zero", "Half", "One"], plot_settings)

    assert _rgba(color_map["Zero"]) == pytest.approx(_rgba(matplotlib.colormaps["plasma"](0.2)))
    assert _rgba(color_map["Half"]) == pytest.approx(_rgba(matplotlib.colormaps["plasma"](0.5)))
    assert _rgba(color_map["One"]) == pytest.approx(_rgba(matplotlib.colormaps["plasma"](0.8)))


def test_semantic_value_colors_override_family_colormap() -> None:
    """Explicit value colors should override family colormap sampling."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "conditions": {"Special": {"family": "category", "value": "x"}},
            "families": {
                "category": {
                    "scale": "ordinal",
                    "colormap": "viridis",
                    "value_colors": {"x": "cyan"},
                }
            },
        }
    )

    assert get_condition_color_map(["Special"], plot_settings)["Special"] == "cyan"


def test_invalid_semantic_condition_color_falls_back(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Invalid semantic condition colors should warn and use missing_color."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "missing_color": "pink",
            "conditions": {"Bad": {"color": "not-a-color"}},
        }
    )

    color_map = get_condition_color_map(["Bad"], plot_settings)

    assert color_map["Bad"] == "pink"
    assert "Invalid condition color for 'Bad'" in caplog.text
    assert "not-a-color" in caplog.text


def test_invalid_semantic_family_colormap_falls_back(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Invalid semantic family colormaps should warn and use missing_color."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "missing_color": "pink",
            "conditions": {"Bad": {"family": "dose", "value": 1.0}},
            "families": {"dose": {"colormap": "not-a-colormap"}},
        }
    )

    color_map = get_condition_color_map(["Bad"], plot_settings)

    assert color_map["Bad"] == "pink"
    assert "Semantic color colormap 'not-a-colormap'" in caplog.text
    assert "Falling back" in caplog.text


def test_semantic_unknown_metadata_uses_missing_color(caplog: pytest.LogCaptureFixture) -> None:
    """Unresolvable semantic metadata should use missing_color with a warning."""
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "missing_color": "pink",
            "conditions": {"Unknown": {"family": "not-configured", "value": 1.0}},
        }
    )

    color_map = get_condition_color_map(["Unknown"], plot_settings)

    assert color_map["Unknown"] == "pink"
    assert "unknown semantic color family" in caplog.text


def test_semantic_no_metadata_uses_default_or_palette() -> None:
    """Labels without metadata should use default_color or palette fallback."""
    labels = ["A", "B"]
    default_settings = PlotSettings(
        semantic_colors={"enabled": True, "default_color": "gray"},
    )
    palette_settings = PlotSettings(semantic_colors={"enabled": True})

    assert get_condition_color_map(labels, default_settings)["A"] == "gray"
    assert get_condition_color_map(labels, palette_settings) == dict(
        zip(labels, get_palette_colors(2, palette_settings))
    )


def test_load_canonical_plot_artifacts_reads_configured_artifacts_only(tmp_path) -> None:
    """Plot artifact loading should ignore non-canonical JSON and extra run directories."""
    analysis_dir = tmp_path / "condition" / "rmsd"
    aggregated_dir = analysis_dir / "aggregated"
    run_1 = analysis_dir / "run_1"
    run_99 = analysis_dir / "run_99"
    run_99.mkdir(parents=True)

    # Replicates are written before the aggregate that summarizes them, the
    # order the pipeline uses; an aggregate older than its replicates is stale.
    ArtifactStore(run_1).write_replicate_result(
        ReplicateArtifact(
            analysis_name="rmsd",
            condition_label="condition",
            replicate=1,
            payload={"metric": 1.1},
        )
    )
    ArtifactStore(aggregated_dir).write_condition_result(
        ConditionArtifact(
            analysis_name="rmsd",
            condition_label="condition",
            replicates=[1],
            payload={"metric": 1.0},
        )
    )
    (run_99 / "result.json").write_text('{"artifact_type": "replicate"}', encoding="utf-8")
    (analysis_dir / "noncanonical_plot.json").write_text('{"metric": 99}', encoding="utf-8")

    loaded = load_canonical_plot_artifacts(analysis_dir, [1])

    assert loaded.condition_artifact is not None
    assert loaded.condition_artifact.payload["metric"] == pytest.approx(1.0)
    assert set(loaded.replicate_artifacts) == {1}
    assert loaded.replicate_artifacts[1].payload["metric"] == pytest.approx(1.1)


def test_load_canonical_plot_artifacts_rejects_condition_model_json(tmp_path) -> None:
    """Non-canonical JSON at the canonical path should fail artifact validation."""
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


def test_load_canonical_plot_artifacts_requires_configured_replicates(tmp_path) -> None:
    """Configured replicate artifacts should be required by default."""
    analysis_dir = tmp_path / "condition" / "rmsd"

    with pytest.raises(ArtifactStoreError, match="Missing canonical replicate artifact"):
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
    values = _finite_numeric_values([1.0, "bad", float("nan"), "2.5", float("inf")])

    np.testing.assert_allclose(values, [1.0, 2.5])


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
