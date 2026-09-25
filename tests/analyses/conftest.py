"""Shared fixtures for the analyses test suite, including the contract harness.

Two kinds of helper live here. The plot audit fixture wraps every plotter so a
figure that draws an uncertainty must also say what it is. The contract harness
gives a plugin author a universe whose answer is known by hand plus a way to run
the real lifecycle without a trajectory on disk, so an end-to-end test is about
five lines.
"""

from __future__ import annotations

import importlib
import types
from pathlib import Path
from typing import Any, Callable, Sequence

import pytest

from polyzymd.analyses.testing import synthetic_universe as make_synthetic_universe

#: Every plugin is on the observable contract, so the framework draws every
#: figure and there is one module left to audit.
_PLOTTER_MODULES = ("polyzymd.analyses.contract_plots",)


def figure_draws_uncertainty(fig: Any) -> bool:
    """Return whether any axes on *fig* draws an error bar or a shaded band."""

    from matplotlib.collections import PolyCollection

    return any(
        any(
            getattr(container, "has_yerr", False) or getattr(container, "has_xerr", False)
            for container in ax.containers
        )
        or any(isinstance(collection, PolyCollection) for collection in ax.collections)
        for ax in fig.axes
    )


def figure_has_uncertainty_footnote(fig: Any) -> bool:
    """Return whether *fig* carries a text naming what its uncertainty is."""

    return any(
        ("95%" in text.get_text() or "SEM" in text.get_text()) and "replicates" in text.get_text()
        for text in fig.texts
    )


@pytest.fixture(autouse=True)
def audit_plot_uncertainty_footnotes(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail any test whose plotter saves an uncertainty figure with no footnote.

    Grossfield et al. (2018) require every figure to describe the meaning and
    basis of its uncertainties. Wrapping ``save_figure`` in every plotter module
    turns the existing body of plot tests into that check on real rendered
    figures. A test that installs its own ``save_figure`` stub bypasses the
    audit; ``tests/analyses/scientific/test_uncertainty_plots.py`` renders those
    paths directly instead.
    """

    pytest.importorskip("matplotlib")

    for module_name in _PLOTTER_MODULES:
        module = importlib.import_module(module_name)
        original = getattr(module, "save_figure", None)
        if original is None:
            continue

        def _checked(
            fig: Any,
            *args: Any,
            _original: Any = original,
            _name: str = module_name,
            **kwargs: Any,
        ) -> Any:
            if figure_draws_uncertainty(fig) and not figure_has_uncertainty_footnote(fig):
                raise AssertionError(
                    f"{_name} saved a figure that draws an uncertainty "
                    "without a footnote saying what it is"
                )
            return _original(fig, *args, **kwargs)

        monkeypatch.setattr(module, "save_figure", _checked)


_DEFAULT_INPUTS = ({"path": "/tmp/topology.pdb", "format": "pdb", "size_bytes": 1, "mtime_ns": 2},)


def make_simulation_config(label: str = "A") -> types.SimpleNamespace:
    """Stand-in carrying only the fields the analysis config hash reads."""
    return types.SimpleNamespace(
        name=label,
        enzyme=types.SimpleNamespace(name="enzyme", pdb_path="/tmp/enzyme.pdb"),
        thermodynamics=types.SimpleNamespace(temperature=300.0, pressure=1.0),
        output=types.SimpleNamespace(
            projects_directory="/tmp/projects",
            effective_scratch_directory="/tmp/scratch",
            naming_template="{name}",
        ),
        substrate=None,
        polymers=None,
    )


@pytest.fixture
def synthetic_universe() -> Any:
    """A five-frame universe whose radius of gyration is exactly 1.0."""
    return make_synthetic_universe()


@pytest.fixture
def serve_replicates(monkeypatch: pytest.MonkeyPatch) -> Callable[..., None]:
    """Serve in-memory universes to the framework in place of trajectories.

    The returned callable takes one universe, or a ``replicate -> universe``
    factory, and optionally ``inputs``: the file identity records the cache
    compares, as a list or a ``replicate -> list`` callable. After it is called,
    every framework entry point reads these universes, with every frame as the
    production window, and no trajectory is opened. It is the only place the
    test suite reaches into how the framework loads a replicate.
    """

    def install(universe: Any, inputs: Any = _DEFAULT_INPUTS) -> None:
        factory = universe if callable(universe) else (lambda replicate: universe)
        current = inputs if callable(inputs) else (lambda replicate: list(inputs))
        _stub_universe_source(monkeypatch, factory, current)

    return install


@pytest.fixture
def run_contract_analysis(
    tmp_path: Path, serve_replicates: Callable[..., None]
) -> Callable[..., Any]:
    """Run a contract analysis for one condition on in-memory universes.

    The returned callable takes the analysis class, its settings, and either a
    universe or a ``replicate -> universe`` factory, and returns the
    ``ConditionArtifact`` that ``orchestrator.run_analysis`` produced.

    Returns
    -------
    Callable[..., ConditionArtifact]
        Harness callable. Keyword arguments: ``label``, ``replicates``,
        ``equilibration``, ``inputs`` (the file identity records the cache
        compares) and ``root``.
    """

    def run(
        analysis_cls: type,
        settings: Any,
        universe: Any,
        *,
        label: str = "A",
        replicates: Sequence[int] = (1, 2, 3),
        equilibration: str = "0ns",
        inputs: Any = _DEFAULT_INPUTS,
        root: Path | None = None,
    ) -> Any:
        from polyzymd.analyses.orchestrator import run_analysis

        base = root or tmp_path
        serve_replicates(universe, inputs)
        return run_analysis(
            analysis_cls(),
            make_condition(label, base, replicates),
            settings,
            equilibration,
            base / "analysis" / label / analysis_cls.name,
        )

    return run


def make_condition(label: str, root: Path, replicates: Sequence[int] = (1, 2, 3)) -> Any:
    """A condition whose simulation config needs no files on disk."""
    from polyzymd.analyses.base import Condition

    return Condition(
        label=label,
        config_path=root / f"{label}.yaml",
        replicates=tuple(replicates),
        sim_config=make_simulation_config(label),
    )


def make_comparison(
    root: Path,
    *,
    labels: Sequence[str] = ("A", "B"),
    replicates: Sequence[int] = (1, 2, 3),
    control: str | None = "A",
    settings: dict[str, Any] | None = None,
    equilibration: str = "0ns",
    **defaults: Any,
) -> Any:
    """A comparison config over conditions whose simulation configs need no files.

    ``settings`` maps analysis name to its settings; ``defaults`` sets
    ``fdr_alpha``, ``ttest_method`` or ``posthoc_method``. Use it together with
    ``serve_replicates``, which makes these conditions load without YAML.
    """
    from types import SimpleNamespace

    from polyzymd.config.comparison import PlotSettings

    plugins = dict(settings or {})
    config = SimpleNamespace(
        name="project",
        source_path=root / "comparison.yaml",
        defaults=SimpleNamespace(equilibration_time=equilibration, **defaults),
        control=control,
        conditions=[
            SimpleNamespace(label=label, config=root / f"{label}.yaml", replicates=list(replicates))
            for label in labels
        ],
        plugins=SimpleNamespace(get=plugins.get, get_enabled_plugins=lambda: list(plugins)),
        plot_settings=PlotSettings(output_dir=root / "figures"),
    )
    config.model_copy = lambda deep=True: config
    return config


def _stub_universe_source(
    monkeypatch: pytest.MonkeyPatch,
    factory: Callable[[int], Any],
    inputs: Callable[[int], list[dict[str, Any]]],
) -> None:
    """Replace the loader, the universe provider and the window the framework uses."""
    from polyzymd.analyses._framework import lifecycle as framework_lifecycle
    from polyzymd.analyses.base import Analysis, Condition
    from polyzymd.analyses.mda import lifecycle
    from polyzymd.analyses.shared.window import TrajectoryWindow

    class _Provider:
        def __init__(self, config: Any, loader: Any = None) -> None:
            self.config = config

        @classmethod
        def from_config(cls, config: Any, loader: Any = None) -> "_Provider":
            return cls(config, loader=loader)

        def load_universe(self, replicate: int) -> Any:
            return factory(replicate)

        def provenance_for(self, replicate: int) -> dict[str, Any]:
            return {"topology": None, "trajectories": inputs(replicate), "warnings": []}

    class _Loader:
        def __init__(self, config: Any) -> None:
            self.config = config

        def get_trajectory_info(self, replicate: int) -> Any:
            files = [Path(entry["path"]) for entry in inputs(replicate) if "path" in entry]
            return types.SimpleNamespace(trajectory_files=files)

    def full_window(self: Any, ctx: Any, replicate: int, loader: Any, universe: Any) -> Any:
        n_frames = len(universe.trajectory)
        return TrajectoryWindow(
            start=0,
            stop=n_frames,
            step=1,
            equilibration_start=0,
            n_frames_total=n_frames,
            n_frames_selected=n_frames,
            timestep_ps=1.0,
            equilibration_ps=0.0,
            equilibration=ctx.equilibration,
        )

    def condition_from_config(cfg: Any) -> Any:
        return Condition(
            cfg.label, Path(cfg.config), tuple(cfg.replicates), make_simulation_config(cfg.label)
        )

    monkeypatch.setattr(Condition, "from_condition_config", staticmethod(condition_from_config))
    monkeypatch.setattr(lifecycle, "UniverseProvider", _Provider)
    monkeypatch.setattr(lifecycle, "build_trajectory_loader", lambda config: _Loader(config))
    monkeypatch.setattr(
        framework_lifecycle, "build_trajectory_loader", lambda config: _Loader(config)
    )
    monkeypatch.setattr(Analysis, "get_trajectory_window", full_window)
