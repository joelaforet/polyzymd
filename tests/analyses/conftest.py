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

_PLOTTER_MODULES = tuple(
    f"polyzymd.analyses.{name}._plotters"
    for name in (
        "rmsd",
        "rmsf",
        "rg",
        "contacts",
        "distances",
        "hydrogen_bonds",
        "secondary_structure",
        "catalytic_triad",
    )
)


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


#: Four atoms on a unit cross. The radius of gyration of this shape is its scale.
CROSS = ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0))

_DEFAULT_INPUTS = ({"path": "/tmp/topology.pdb", "format": "pdb", "size_bytes": 1, "mtime_ns": 2},)


def make_synthetic_universe(scale: float = 1.0, n_frames: int = 5) -> Any:
    """Build an in-memory universe of four unit-mass atoms on a cross.

    Parameters
    ----------
    scale : float, optional
        Distance of each atom from the origin, by default 1.0. The radius of
        gyration of the group is exactly this value.
    n_frames : int, optional
        Number of identical frames, by default 5.

    Returns
    -------
    MDAnalysis.Universe
        Universe backed by ``MemoryReader``.
    """
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("masses", [1.0] * 4)
    positions = np.asarray(CROSS, dtype=np.float32) * scale
    universe.load_new(np.stack([positions] * n_frames), format=MemoryReader)
    return universe


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
def run_contract_analysis(tmp_path: Path) -> Callable[..., Any]:
    """Run a contract analysis end to end on in-memory universes.

    The returned callable takes the generated analysis class, its settings, and
    either a universe or a ``replicate -> universe`` factory, and returns the
    ``ConditionArtifact`` that ``AnalysisLifecycle.run_analysis`` produced. The
    trajectory loader and universe provider are replaced, so no files are read.

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
        inputs: Sequence[dict[str, Any]] = _DEFAULT_INPUTS,
        root: Path | None = None,
    ) -> Any:
        from polyzymd.analyses._framework.contexts import Condition
        from polyzymd.analyses._framework.lifecycle import AnalysisLifecycle

        base = root or tmp_path
        factory = universe if callable(universe) else (lambda replicate: universe)
        condition = Condition(
            label=label,
            config_path=base / f"{label}.yaml",
            replicates=tuple(replicates),
            sim_config=make_simulation_config(label),
        )
        stub_cls = _stubbed(analysis_cls, factory, tuple(inputs))
        return AnalysisLifecycle(stub_cls()).run_analysis(
            condition,
            settings,
            equilibration,
            base / "analysis" / label / stub_cls.name,
        )

    return run


def _stubbed(analysis_cls: type, factory: Callable[[int], Any], inputs: tuple) -> type:
    """Subclass an analysis with the loader and universe provider replaced."""
    from polyzymd.analyses.shared.window import TrajectoryWindow

    class _Provider:
        def __init__(self, config: Any, loader: Any = None) -> None:
            self.config = config

        def load_universe(self, replicate: int) -> Any:
            return factory(replicate)

        def provenance_for(self, replicate: int) -> dict[str, Any]:
            return {"topology": None, "trajectories": list(inputs), "warnings": []}

    class _Loader:
        def __init__(self, config: Any) -> None:
            self.config = config

    class _Stubbed(analysis_cls):  # type: ignore[valid-type, misc]
        def _trajectory_loader_factory(self) -> type:
            return _Loader

        def _mda_universe_provider_factory(self) -> type:
            return _Provider

        def get_trajectory_window(
            self, ctx: Any, replicate: int, loader: Any, universe: Any
        ) -> TrajectoryWindow:
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

    return _Stubbed
