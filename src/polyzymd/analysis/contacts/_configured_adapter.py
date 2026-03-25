"""Config-driven adapter for contact analysis.

This adapter wraps the universe-driven ContactAnalyzer/ParallelContactAnalyzer
behind the BaseAnalyzer interface, enabling registry-driven creation.

The original ContactAnalyzer is universe-driven (run(universe, ...)) which is a
different abstraction from config-driven analyzers. This adapter bridges the
gap without modifying the original ContactAnalyzer.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Sequence

from polyzymd.analysis.core.constants import DEFAULT_CONTACT_CUTOFF
from polyzymd.analysis.core.registry import AnalyzerRegistry, BaseAnalysisSettings, BaseAnalyzer

if TYPE_CHECKING:
    from polyzymd.analysis.results.base import BaseAnalysisResult
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger(__name__)


@AnalyzerRegistry.register("contacts")
class ConfiguredContactsAnalyzer(BaseAnalyzer):
    """Config-driven adapter for contact analysis.

    This wraps ParallelContactAnalyzer behind the BaseAnalyzer interface.
    The original ContactAnalyzer remains unchanged for direct use.

    Notes
    -----
    This adapter intentionally maps only ``polymer_selection``,
    ``protein_selection``, and ``cutoff`` from analysis settings. Binding
    preference and grouping behavior are configured by CLI orchestration and
    comparator layers, not by this minimal wrapper.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration for trajectory loading.
    polymer_selection : str, optional
        MDAnalysis selection for polymer atoms.
    protein_selection : str, optional
        MDAnalysis selection for protein atoms.
    cutoff : float, optional
        Contact distance cutoff in Angstroms.
    equilibration : str, optional
        Equilibration time to skip.
    """

    def __init__(
        self,
        sim_config: "SimulationConfig",
        polymer_selection: str = "chainID C",
        protein_selection: str = "protein",
        cutoff: float = DEFAULT_CONTACT_CUTOFF,
        equilibration: str = "0ns",
    ) -> None:
        self._sim_config = sim_config
        self._polymer_selection = polymer_selection
        self._protein_selection = protein_selection
        self._cutoff = cutoff
        self._equilibration = equilibration

    @classmethod
    def analysis_type(cls) -> str:
        """Return the unique identifier for this analyzer.

        Returns
        -------
        str
            Analysis type identifier.
        """
        return "contacts"

    @classmethod
    def from_config(
        cls,
        analysis_settings: BaseAnalysisSettings,
        sim_config: "SimulationConfig",
        equilibration: str = "0ns",
    ) -> "ConfiguredContactsAnalyzer":
        """Create a configured contacts analyzer from analysis settings.

        Parameters
        ----------
        analysis_settings : BaseAnalysisSettings
            Contacts-compatible settings object.
        sim_config : SimulationConfig
            Simulation configuration for trajectory loading.
        equilibration : str, optional
            Equilibration time to skip.

        Returns
        -------
        ConfiguredContactsAnalyzer
            Configured contacts analyzer.
        """
        return cls(
            sim_config=sim_config,
            polymer_selection=getattr(analysis_settings, "polymer_selection", "chainID C"),
            protein_selection=getattr(analysis_settings, "protein_selection", "protein"),
            cutoff=getattr(analysis_settings, "cutoff", DEFAULT_CONTACT_CUTOFF),
            equilibration=equilibration,
        )

    def compute(self, replicate: int, **kwargs) -> "BaseAnalysisResult":
        """Run contact analysis for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number to analyze.
        **kwargs
            Additional analyzer-specific parameters.

        Returns
        -------
        BaseAnalysisResult
            Contact analysis result for the replicate.
        """
        from polyzymd.analysis.common.selectors import MDAnalysisSelector
        from polyzymd.analysis.contacts import ParallelContactAnalyzer
        from polyzymd.analysis.core.loader import (
            TrajectoryLoader,
            convert_time,
            parse_time_string,
            time_to_frame,
        )

        loader = TrajectoryLoader(self._sim_config)
        universe = loader.load_universe(replicate)

        eq_value, eq_unit = parse_time_string(self._equilibration)
        timestep = loader.get_timestep(replicate)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")
        start_frame = time_to_frame(eq_time_ps, "ps", timestep, "ps")

        target_selector = MDAnalysisSelector(self._protein_selection)
        query_selector = MDAnalysisSelector(self._polymer_selection)

        analyzer = ParallelContactAnalyzer(
            target_selector=target_selector,
            query_selector=query_selector,
            cutoff=self._cutoff,
        )

        return analyzer.run(universe, start=start_frame)

    def compute_aggregated(self, replicates: Sequence[int], **kwargs) -> "BaseAnalysisResult":
        """Run aggregated contact analysis across replicates.

        Parameters
        ----------
        replicates : Sequence[int]
            List of replicate numbers to analyze.
        **kwargs
            Additional analyzer-specific parameters.

        Returns
        -------
        BaseAnalysisResult
            Aggregated contact analysis result.
        """
        from polyzymd.analysis.contacts.aggregator import aggregate_contact_results

        results = [self.compute(rep, **kwargs) for rep in replicates]
        return aggregate_contact_results(results)
