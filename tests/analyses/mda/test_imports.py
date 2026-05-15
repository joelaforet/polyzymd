"""Tests for MDAnalysis extension-layer public imports."""

from __future__ import annotations

import importlib
import sys
from typing import Any


def test_public_facade_reexports_primitives() -> None:
    """The package facade should expose the stable extension-layer primitives."""

    from polyzymd.analyses import mda
    from polyzymd.analyses.mda.base import (
        MDA_EXTENSION_API_VERSION,
        AnalysisBaseLike,
        MDAnalysisExtensionError,
        MDARunKwargs,
    )

    assert mda.MDA_EXTENSION_API_VERSION == MDA_EXTENSION_API_VERSION == "1"
    assert mda.AnalysisBaseLike is AnalysisBaseLike
    assert mda.MDAnalysisExtensionError is MDAnalysisExtensionError
    assert mda.MDARunKwargs is MDARunKwargs
    assert set(mda.__all__) == {
        "MDA_EXTENSION_API_VERSION",
        "AnalysisBaseLike",
        "MDAnalysisExtensionError",
        "MDARunKwargs",
    }


def test_import_does_not_load_heavy_simulation_modules() -> None:
    """Importing the extension primitives should not import MDAnalysis or OpenMM."""

    heavy_modules = ("MDAnalysis", "openmm", "openff", "parmed", "pdbfixer")
    initially_loaded = {name for name in heavy_modules if name in sys.modules}

    importlib.import_module("polyzymd.analyses.mda")
    importlib.import_module("polyzymd.analyses.mda.base")

    for module_name in heavy_modules:
        if module_name not in initially_loaded:
            assert module_name not in sys.modules


def test_analysis_base_like_protocol_accepts_run_results_object() -> None:
    """The protocol should match objects with ``results`` and chainable ``run``."""

    from polyzymd.analyses.mda import AnalysisBaseLike, MDARunKwargs

    class FakeAnalysisBase:
        """Minimal structural stand-in for an MDAnalysis analysis object."""

        results: dict[str, Any]

        def __init__(self) -> None:
            self.results = {}
            self.run_kwargs: dict[str, Any] = {}

        def run(self, **kwargs: Any) -> "FakeAnalysisBase":
            """Store run keyword arguments and return ``self``."""

            self.run_kwargs = dict(kwargs)
            return self

    kwargs = MDARunKwargs(start=1, stop=10, step=2, verbose=False)
    analysis = FakeAnalysisBase()

    assert isinstance(analysis, AnalysisBaseLike)
    assert analysis.run(**kwargs) is analysis
    assert analysis.run_kwargs == kwargs


def test_extension_error_is_runtime_error() -> None:
    """Extension failures should use a RuntimeError subclass."""

    from polyzymd.analyses.mda import MDAnalysisExtensionError

    assert issubclass(MDAnalysisExtensionError, RuntimeError)
