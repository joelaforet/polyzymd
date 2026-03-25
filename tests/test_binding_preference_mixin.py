"""Regression tests for BindingPreferenceComparatorMixin.

Tests verify the shared load-or-compute workflow, hook dispatch, and
settings resolution used by both BindingFreeEnergyComparator and
PolymerAffinityScoreComparator.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import patch

from polyzymd.compare.comparators.binding_free_energy import BindingFreeEnergyComparator
from polyzymd.compare.comparators.polymer_affinity import PolymerAffinityScoreComparator
from polyzymd.compare.settings import (
    BindingFreeEnergyAnalysisSettings,
    PolymerAffinityScoreSettings,
)


def _make_condition(config_path: Path | None = None) -> SimpleNamespace:
    """Build a minimal condition-like object for comparator hooks.

    Parameters
    ----------
    config_path : Path | None, optional
        Path for the condition config field, by default None.

    Returns
    -------
    SimpleNamespace
        Object with label, config, and replicates attributes.
    """
    return SimpleNamespace(
        label="TestCond",
        config=config_path or Path("/fake/config.yaml"),
        replicates=[1, 2, 3],
    )


def _make_config(cond: Any, contacts_settings: Any | None = None) -> SimpleNamespace:
    """Build a minimal comparison-like config object.

    Parameters
    ----------
    cond : Any
        Condition object to include.
    contacts_settings : Any | None, optional
        Optional contacts settings for fallback tests, by default None.

    Returns
    -------
    SimpleNamespace
        Object with comparator-required config fields.
    """
    analysis_settings = {}
    if contacts_settings is not None:
        analysis_settings["contacts"] = contacts_settings

    return SimpleNamespace(
        source_path=None,
        name="test",
        conditions=[cond],
        control=None,
        analysis_settings=analysis_settings,
    )


def _make_sim_config(temperature: float = 300.0) -> SimpleNamespace:
    """Build a minimal simulation config object for lazy loading.

    Parameters
    ----------
    temperature : float, optional
        Thermodynamic temperature in Kelvin, by default 300.0.

    Returns
    -------
    SimpleNamespace
        Simulation-like object with thermodynamics.temperature.
    """
    return SimpleNamespace(thermodynamics=SimpleNamespace(temperature=temperature))


class TestLoadOrComputeCacheHit:
    """Regression tests for cache-hit load-or-compute behavior.

    Parameters
    ----------
    None
        This class validates returned payload shape and cache short-circuiting.
    """

    def test_bfe_cache_hit_returns_base_keys(self, tmp_path: Path) -> None:
        """Return core keys for BFE when cached data exists.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        bp_result = object()
        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=bp_result,
            ),
        ):
            result = comparator._load_or_compute(cond, recompute=False)

        assert set(result.keys()) == {"bp_result", "temperature_K", "cond_label", "config_path"}
        assert result["bp_result"] is bp_result
        assert result["temperature_K"] == 300.0

    def test_pa_cache_hit_returns_extra_keys(self, tmp_path: Path) -> None:
        """Return extra PA keys in addition to base keys on cache hit.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = PolymerAffinityScoreComparator(
            config,
            PolymerAffinityScoreSettings(),
            equilibration="10ns",
        )

        bp_result = object()
        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=bp_result,
            ),
        ):
            result = comparator._load_or_compute(cond, recompute=False)

        assert set(result.keys()) == {
            "bp_result",
            "temperature_K",
            "cond_label",
            "config_path",
            "analysis_dir",
            "cond",
        }
        assert result["analysis_dir"] == analysis_dir
        assert result["cond"] is cond

    def test_cache_hit_skips_compute(self, tmp_path: Path) -> None:
        """Skip on-demand compute call when cache returns a result.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=object(),
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".compute_condition_binding_preference"
            ) as compute_mock,
        ):
            comparator._load_or_compute(cond, recompute=False)

        compute_mock.assert_not_called()


class TestLoadOrComputeCacheMiss:
    """Regression tests for cache-miss load-or-compute behavior.

    Parameters
    ----------
    None
        This class validates compute fallback and failure pathways.
    """

    def test_bfe_compute_on_cache_miss(self, tmp_path: Path) -> None:
        """Compute and return BFE binding preference when cache misses.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("MODEL\nENDMDL\n")
        computed_bp = object()

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=None,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".compute_condition_binding_preference",
                return_value=computed_bp,
            ),
        ):
            result = comparator._load_or_compute(cond, recompute=False)

        assert result["bp_result"] is computed_bp

    def test_pa_compute_on_cache_miss_has_extra_keys(self, tmp_path: Path) -> None:
        """Include PA extra keys when computed after cache miss.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = PolymerAffinityScoreComparator(
            config,
            PolymerAffinityScoreSettings(),
            equilibration="10ns",
        )

        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("MODEL\nENDMDL\n")
        computed_bp = object()

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=None,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".compute_condition_binding_preference",
                return_value=computed_bp,
            ),
        ):
            result = comparator._load_or_compute(cond, recompute=False)

        assert result["bp_result"] is computed_bp
        assert result["analysis_dir"] == analysis_dir
        assert result["cond"] is cond

    def test_missing_enzyme_pdb_returns_none_bp_result(self, tmp_path: Path) -> None:
        """Return None bp_result when enzyme PDB cannot be resolved.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=None,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers.resolve_enzyme_pdb",
                return_value=None,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".compute_condition_binding_preference"
            ) as compute_mock,
        ):
            result = comparator._load_or_compute(cond, recompute=False)

        assert result["bp_result"] is None
        compute_mock.assert_not_called()

    def test_compute_disabled_returns_none_bp_result(self, tmp_path: Path) -> None:
        """Return None bp_result when compute is disabled in settings.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        analysis_settings = BindingFreeEnergyAnalysisSettings(compute_binding_preference=False)
        comparator = BindingFreeEnergyComparator(
            config,
            analysis_settings,
            equilibration="10ns",
        )

        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference",
                return_value=None,
            ),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".compute_condition_binding_preference"
            ) as compute_mock,
        ):
            result = comparator._load_or_compute(cond, recompute=False)

        assert result["bp_result"] is None
        compute_mock.assert_not_called()


class TestLoadOrComputeRecompute:
    """Regression tests for recompute control flow.

    Parameters
    ----------
    None
        This class validates that recompute bypasses cache loading.
    """

    def test_recompute_true_skips_cache(self, tmp_path: Path) -> None:
        """Skip cache loader when recompute flag is True.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        analysis_settings = BindingFreeEnergyAnalysisSettings(compute_binding_preference=False)
        comparator = BindingFreeEnergyComparator(
            config,
            analysis_settings,
            equilibration="10ns",
        )

        sim_config = _make_sim_config(300.0)
        analysis_dir = tmp_path / "contacts"

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch.object(comparator, "_resolve_condition_output_dir", return_value=analysis_dir),
            patch(
                "polyzymd.compare.comparators._binding_preference_helpers"
                ".try_load_cached_binding_preference"
            ) as cached_mock,
        ):
            comparator._load_or_compute(cond, recompute=True)

        cached_mock.assert_not_called()


class TestResolveComputeSettings:
    """Regression tests for compute-settings resolution and fallback.

    Parameters
    ----------
    None
        This class validates direct, fallback, and default setting resolution.
    """

    def test_settings_from_analysis_settings(self) -> None:
        """Prefer attributes set on comparator analysis settings.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        analysis_settings = BindingFreeEnergyAnalysisSettings(
            enzyme_pdb_for_sasa="enzyme_from_bfe.pdb"
        )
        comparator = BindingFreeEnergyComparator(
            config,
            analysis_settings,
            equilibration="10ns",
        )

        settings = comparator._resolve_compute_settings()

        assert settings["enzyme_pdb_for_sasa"] == "enzyme_from_bfe.pdb"

    def test_fallback_to_contacts_settings(self) -> None:
        """Fallback to contacts settings when analysis value is None.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        contacts_settings = SimpleNamespace(enzyme_pdb_for_sasa="enzyme_from_contacts.pdb")
        config = _make_config(cond, contacts_settings=contacts_settings)
        analysis_settings = BindingFreeEnergyAnalysisSettings(enzyme_pdb_for_sasa=None)
        comparator = BindingFreeEnergyComparator(
            config,
            analysis_settings,
            equilibration="10ns",
        )

        settings = comparator._resolve_compute_settings()

        assert settings["enzyme_pdb_for_sasa"] == "enzyme_from_contacts.pdb"

    def test_defaults_when_no_settings(self) -> None:
        """Use hard-coded defaults when no settings provide values.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        analysis_settings = SimpleNamespace(
            enzyme_pdb_for_sasa=None,
            surface_exposure_threshold=None,
            include_default_aa_groups=None,
            protein_groups=None,
            protein_partitions=None,
            polymer_type_selections=None,
        )
        comparator = BindingFreeEnergyComparator(
            config,
            analysis_settings,
            equilibration="10ns",
        )

        settings = comparator._resolve_compute_settings()

        assert settings["surface_exposure_threshold"] == 0.2
        assert settings["include_default_aa_groups"] is True


class TestFindContactsAnalysisDir:
    """Regression tests for contacts analysis directory resolution.

    Parameters
    ----------
    None
        This class validates direct-return and delegated discovery behavior.
    """

    def test_returns_condition_output_dir_when_provided(self, tmp_path: Path) -> None:
        """Return condition output directory when provided.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory provided by pytest.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        condition_output_dir = tmp_path / "analysis" / "contacts"
        result = comparator._find_contacts_analysis_dir(
            _make_sim_config(),
            cond,
            condition_output_dir=condition_output_dir,
        )

        assert result == condition_output_dir

    def test_delegates_to_find_analysis_dir_when_no_condition_dir(self) -> None:
        """Delegate lookup to find_analysis_dir when no condition dir exists.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        sim_config = _make_sim_config()
        expected_dir = Path("/resolved/analysis/contacts")

        with patch(
            "polyzymd.compare.comparators._utils.find_analysis_dir",
            return_value=expected_dir,
        ) as find_mock:
            result = comparator._find_contacts_analysis_dir(
                sim_config,
                cond,
                condition_output_dir=None,
            )

        assert result == expected_dir
        find_mock.assert_called_once_with(
            sim_config,
            analysis_subdir="analysis/contacts",
            cond_config_path=Path(cond.config),
        )


class TestHookDispatch:
    """Regression tests for subclass hook dispatch behavior.

    Parameters
    ----------
    None
        This class validates the two subclass-overridden hook methods.
    """

    def test_bfe_settings_label(self) -> None:
        """Return binding_free_energy settings label in BFE comparator.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        assert comparator._binding_preference_settings_label() == "binding_free_energy"

    def test_pa_settings_label(self) -> None:
        """Return polymer_affinity settings label in PA comparator.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = PolymerAffinityScoreComparator(
            config,
            PolymerAffinityScoreSettings(),
            equilibration="10ns",
        )

        assert comparator._binding_preference_settings_label() == "polymer_affinity"

    def test_bfe_extra_data_empty(self) -> None:
        """Return empty extra data for BFE condition payloads.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = BindingFreeEnergyComparator(
            config,
            BindingFreeEnergyAnalysisSettings(),
            equilibration="10ns",
        )

        extra = comparator._binding_preference_extra_condition_data(cond, Path("/tmp"), object())

        assert extra == {}

    def test_pa_extra_data_has_analysis_dir_and_cond(self) -> None:
        """Return analysis_dir and cond extra keys for PA payloads.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        cond = _make_condition()
        config = _make_config(cond)
        comparator = PolymerAffinityScoreComparator(
            config,
            PolymerAffinityScoreSettings(),
            equilibration="10ns",
        )

        analysis_dir = Path("/tmp/contacts")
        extra = comparator._binding_preference_extra_condition_data(cond, analysis_dir, object())

        assert extra == {"analysis_dir": analysis_dir, "cond": cond}
