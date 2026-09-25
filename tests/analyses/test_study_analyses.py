"""Analyses that live in a study folder rather than in the PolyzyMD source tree.

A study is a directory marked by ``study.yaml``. These tests pin how its
``analyses/`` folder is found, what gets registered from it, which names it may
not use, and the public helpers a study author relies on: ``analyze()`` taking a
plugin object, ``load_replicate`` and ``polyzymd.analyses.testing``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from textwrap import dedent
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses import discovery
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.config.study import analyses_directory, find_study_root

PLUGIN_SOURCE = """
from typing import ClassVar

from pydantic import BaseModel

from polyzymd.analyses import Observable, iter_frames


class EndToEndSettings(BaseModel):
    selection: str = "all"


class EndToEnd:
    name: ClassVar[str] = "{name}"
    Settings: ClassVar[type[BaseModel]] = EndToEndSettings
    references: ClassVar[tuple[str, ...]] = ()

    def compute(self, universe, frames, settings):
        group = universe.select_atoms(settings.selection)
        values = [float(group.radius_of_gyration()) for _ in iter_frames(universe, frames)]
        return [Observable(name="end_to_end", kind="mean_of_timeseries", unit="A", values=values)]
"""


@pytest.fixture(autouse=True)
def _forget_registered_analyses():
    """Leave discovery and ``sys.modules`` as each test found them."""
    discovery.clear_cache()
    yield
    discovery.clear_cache()
    for name in [module for module in sys.modules if module.startswith("polyzymd_study_")]:
        del sys.modules[name]


def _study(tmp_path: Path, **files: str) -> Path:
    """Write a study with the given files in its analyses folder."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "study.yaml").write_text("name: probe\n")
    analyses = tmp_path / "analyses"
    analyses.mkdir()
    for filename, source in files.items():
        path = analyses / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(dedent(source))
    return tmp_path


class TestFindingTheStudy:
    """The study root and its analyses folder are found from anywhere inside it."""

    def test_root_is_found_from_a_nested_file(self, tmp_path: Path) -> None:
        root = _study(tmp_path)
        nested = root / "comparisons" / "c1" / "comparison.yaml"
        nested.parent.mkdir(parents=True)
        nested.write_text("name: c1\n")

        assert find_study_root(nested) == root.resolve()
        assert analyses_directory(nested) == (root / "analyses").resolve()

    def test_study_yaml_can_name_another_folder(self, tmp_path: Path) -> None:
        (tmp_path / "study.yaml").write_text("name: probe\nanalyses: code/analyses\n")
        (tmp_path / "code" / "analyses").mkdir(parents=True)

        assert analyses_directory(tmp_path) == (tmp_path / "code" / "analyses").resolve()

    def test_without_a_study_the_folder_beside_the_file_is_used(self, tmp_path: Path) -> None:
        (tmp_path / "analyses").mkdir()
        comparison = tmp_path / "comparison.yaml"
        comparison.write_text("name: c\n")

        assert find_study_root(comparison) is None
        assert analyses_directory(comparison) == (tmp_path / "analyses").resolve()

    def test_no_folder_means_no_study_analyses(self, tmp_path: Path) -> None:
        assert analyses_directory(tmp_path / "comparison.yaml") is None

    def test_unknown_study_keys_are_rejected(self, tmp_path: Path) -> None:
        from polyzymd.config.study import StudyConfig

        (tmp_path / "study.yaml").write_text("name: probe\nanalysis: typo\n")

        with pytest.raises(ValueError, match="analysis"):
            StudyConfig.from_yaml(tmp_path / "study.yaml")


class TestLoadingTheFolder:
    """Every analysis module in the folder is registered, and nothing else."""

    def test_a_bare_plugin_class_is_registered_without_contract_analysis(
        self, tmp_path: Path
    ) -> None:
        root = _study(tmp_path, **{"end_to_end.py": PLUGIN_SOURCE.format(name="end_to_end")})

        assert discovery.load_analysis_directory(root / "analyses") == ["end_to_end"]
        analysis = discovery.get_analysis("end_to_end")
        assert type(analysis.plugin).__name__ == "EndToEnd"
        assert "end_to_end" in discovery.list_all_names()

    def test_a_wrapped_plugin_is_registered_once(self, tmp_path: Path) -> None:
        source = PLUGIN_SOURCE.format(name="wrapped") + (
            "\nfrom polyzymd.analyses import contract_analysis\n"
            "EndToEndAnalysis = contract_analysis(EndToEnd)\n"
        )
        root = _study(tmp_path, **{"wrapped.py": source})

        assert discovery.load_analysis_directory(root / "analyses") == ["wrapped"]

    def test_helpers_and_tests_are_not_imported_as_analyses(self, tmp_path: Path) -> None:
        plugin = PLUGIN_SOURCE.format(name="uses_helper").replace(
            "from polyzymd.analyses import Observable, iter_frames",
            "from polyzymd.analyses import Observable, iter_frames\nfrom ._units import ANGSTROM",
        )
        root = _study(
            tmp_path,
            **{
                "uses_helper.py": plugin,
                "_units.py": "ANGSTROM = 'A'\n",
                "test_uses_helper.py": "raise RuntimeError('a test module was imported')\n",
                "tests/test_more.py": "raise RuntimeError('a tests folder was imported')\n",
            },
        )

        assert discovery.load_analysis_directory(root / "analyses") == ["uses_helper"]

    def test_an_imported_builtin_is_not_registered_again(self, tmp_path: Path) -> None:
        source = "from polyzymd.analyses.rg import Rg, RgAnalysis  # noqa: F401\n"
        root = _study(tmp_path, **{"reuses_rg.py": source})

        assert discovery.load_analysis_directory(root / "analyses") == []

    def test_a_builtin_name_is_refused(self, tmp_path: Path) -> None:
        root = _study(tmp_path, **{"my_rg.py": PLUGIN_SOURCE.format(name="rg")})

        with pytest.raises(PluginContractError, match="built-in PolyzyMD analysis"):
            discovery.load_analysis_directory(root / "analyses")

    def test_two_files_with_one_name_are_refused(self, tmp_path: Path) -> None:
        root = _study(
            tmp_path,
            **{
                "a.py": PLUGIN_SOURCE.format(name="same"),
                "b.py": PLUGIN_SOURCE.format(name="same"),
            },
        )

        with pytest.raises(PluginContractError, match="collision"):
            discovery.load_analysis_directory(root / "analyses")

    def test_an_import_error_names_the_file(self, tmp_path: Path) -> None:
        root = _study(tmp_path, **{"broken.py": "import not_a_real_module_xyz\n"})

        with pytest.raises(PluginContractError, match="broken.py"):
            discovery.load_analysis_directory(root / "analyses")

    def test_two_studies_in_one_process_do_not_share_modules(self, tmp_path: Path) -> None:
        first = _study(tmp_path / "one", **{"first.py": PLUGIN_SOURCE.format(name="first")})
        second = _study(tmp_path / "two", **{"first.py": PLUGIN_SOURCE.format(name="second")})

        discovery.load_analysis_directory(first / "analyses")
        discovery.load_analysis_directory(second / "analyses")

        assert {"first", "second"} <= set(discovery.list_all_names())


class TestComparisonConfig:
    """A comparison inside a study validates the study's analyses by name."""

    def _comparison(self, root: Path, settings: str) -> Path:
        path = root / "comparisons" / "c1" / "comparison.yaml"
        path.parent.mkdir(parents=True)
        path.write_text(
            "name: c1\n"
            "conditions:\n"
            "  - label: A\n"
            "    config: ../../conditions/a/config.yaml\n"
            "    replicates: [1]\n"
            f"plugins:\n  end_to_end:\n{settings}"
        )
        return path

    def test_study_settings_are_validated(self, tmp_path: Path) -> None:
        from polyzymd.config.comparison import ComparisonConfig

        root = _study(tmp_path, **{"end_to_end.py": PLUGIN_SOURCE.format(name="end_to_end")})
        config = ComparisonConfig.from_yaml(self._comparison(root, "    selection: protein\n"))

        assert config.plugins.get("end_to_end").selection == "protein"

    def test_a_bad_setting_is_an_error(self, tmp_path: Path) -> None:
        from polyzymd.config.comparison import ComparisonConfig

        root = _study(tmp_path, **{"end_to_end.py": PLUGIN_SOURCE.format(name="end_to_end")})
        path = self._comparison(root, "    selection: [not, a, string]\n")

        with pytest.raises(ValueError, match="selection"):
            ComparisonConfig.from_yaml(path)


class TestPluginObjects:
    """A plugin defined in a script is accepted wherever a name is."""

    def test_analyze_resolves_a_plugin_object(self) -> None:
        from polyzymd.analyses.protocols import _analysis_class

        namespace: dict = {}
        exec(PLUGIN_SOURCE.format(name="scripted"), namespace)
        analysis = _analysis_class(namespace["EndToEnd"])

        assert analysis.name == "scripted"
        assert discovery.get_analysis("scripted") is analysis

    def test_a_broken_plugin_object_is_a_protocol_error(self) -> None:
        from polyzymd.analyses.exceptions import ProtocolError
        from polyzymd.analyses.protocols import _analysis_class

        class NoCompute:
            name = "no_compute"

        with pytest.raises(ProtocolError, match="AnalysisProtocol"):
            _analysis_class(NoCompute)


class TestInMemoryTesting:
    """``run_in_memory`` gives the numbers the framework would report."""

    def test_known_answer_over_three_replicates(self) -> None:
        from polyzymd.analyses.testing import run_in_memory, synthetic_universe

        namespace: dict = {}
        exec(PLUGIN_SOURCE.format(name="in_memory"), namespace)
        aggregates = run_in_memory(namespace["EndToEnd"], {}, synthetic_universe(scale=2.0))

        assert aggregates[0].n_replicates == 3
        assert aggregates[0].mean == pytest.approx(2.0)
        assert aggregates[0].sem == pytest.approx(0.0)

    def test_distinct_replicates_give_a_spread(self) -> None:
        from polyzymd.analyses.testing import run_in_memory, synthetic_universe

        namespace: dict = {}
        exec(PLUGIN_SOURCE.format(name="spread"), namespace)
        universes = [synthetic_universe(scale=scale) for scale in (1.0, 2.0, 3.0)]
        aggregates = run_in_memory(namespace["EndToEnd"], {}, universes)

        assert aggregates[0].replicate_values == pytest.approx([1.0, 2.0, 3.0])
        assert aggregates[0].sem == pytest.approx(1.0 / 3.0**0.5)


class TestLoadReplicate:
    """``load_replicate`` returns the universe and window a plugin receives."""

    def test_segments_are_joined_and_equilibration_discarded(self, tmp_path: Path) -> None:
        import MDAnalysis as mda
        import numpy as np

        from polyzymd.analyses import load_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        topology = mda.Universe.empty(3, trajectory=True)
        topology.add_TopologyAttr("names", ["C1", "C2", "C3"])
        topology.add_TopologyAttr("resnames", ["MOL"])
        topology.atoms.positions = np.zeros((3, 3))
        topology.atoms.write(str(run_dir / "solvated_system.pdb"))
        for segment in range(2):
            directory = run_dir / f"production_{segment}"
            directory.mkdir()
            with mda.Writer(
                str(directory / f"production_{segment}_trajectory.dcd"),
                n_atoms=3,
                dt=1.0,
                istart=5 * segment,
                nsavc=1,
            ) as writer:
                for _ in range(5):
                    writer.write(topology.atoms)

        config = MagicMock()
        config.engine = "openmm"
        config.get_working_directory.side_effect = lambda replicate: tmp_path / f"run_{replicate}"
        config.output.effective_scratch_directory = tmp_path

        universe, frames = load_replicate(config, 1, equilibration="3ps")

        assert len(universe.trajectory) == 10
        assert (frames.start, frames.stop, frames.step) == (3, 10, 1)
        assert frames.run_kwargs()["start"] == 3
