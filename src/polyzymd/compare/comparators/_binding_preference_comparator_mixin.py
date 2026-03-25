"""Shared load-or-compute helpers for binding preference comparators.

This mixin centralizes the common binding preference data flow used by
comparators that depend on cached or on-demand binding preference results.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from polyzymd.compare.config import ConditionConfig

logger = logging.getLogger("polyzymd.compare")


class BindingPreferenceComparatorMixin:
    """Shared behavior for binding preference-backed comparators.

    Provides a unified load-or-compute workflow with extension hooks for
    subclass-specific metadata and warning-label customization.
    """

    def _binding_preference_settings_label(self) -> str:
        """Return the settings section name for warning messages.

        Returns
        -------
        str
            Settings section label used in warning text.
        """
        return "binding_preference"

    def _binding_preference_extra_condition_data(
        self,
        cond: "ConditionConfig",
        analysis_dir: Path | None,
        bp_result: Any,
    ) -> dict[str, Any]:
        """Return extra keys for condition data dictionaries.

        Parameters
        ----------
        cond : ConditionConfig
            Condition being processed.
        analysis_dir : Path | None
            Resolved contacts analysis directory.
        bp_result : Any
            Loaded or computed binding preference result.

        Returns
        -------
        dict[str, Any]
            Extra key-value pairs to merge into the condition data.
        """
        return {}

    def _load_or_compute(self, cond: "ConditionConfig", recompute: bool) -> dict[str, Any]:
        """Load cached binding preference data or compute it on demand.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to load or compute.
        recompute : bool
            If True, skip cache and recompute.

        Returns
        -------
        dict[str, Any]
            Condition data including binding preference result and metadata.
        """
        from polyzymd.compare.comparators._binding_preference_helpers import (
            compute_condition_binding_preference,
            resolve_enzyme_pdb,
        )
        from polyzymd.config.schema import SimulationConfig

        settings_label = self._binding_preference_settings_label()
        logger.info(f"Loading binding preference for: {cond.label}")

        sim_config = SimulationConfig.from_yaml(cond.config)
        temperature_K = float(sim_config.thermodynamics.temperature)

        condition_output_dir = self._resolve_condition_output_dir(cond.label, "contacts")
        analysis_dir = self._find_contacts_analysis_dir(
            sim_config, cond, condition_output_dir=condition_output_dir
        )

        bp_result = None
        if not recompute:
            bp_result = self._try_load_cached_binding_preference(cond, analysis_dir)

        if bp_result is not None:
            logger.info(f"  Loaded binding preference for {cond.label} at {temperature_K} K")
            base_data = {
                "bp_result": bp_result,
                "temperature_K": temperature_K,
                "cond_label": cond.label,
                "config_path": str(cond.config),
            }
            base_data.update(
                self._binding_preference_extra_condition_data(cond, analysis_dir, bp_result)
            )
            return base_data

        compute_enabled = getattr(self.analysis_settings, "compute_binding_preference", True)
        if compute_enabled:
            logger.info(f"  No cached data for {cond.label}, computing binding preference...")
            settings = self._resolve_compute_settings()

            enzyme_pdb = resolve_enzyme_pdb(
                enzyme_pdb_setting=settings["enzyme_pdb_for_sasa"],
                source_path=self.config.source_path,
                sim_config=sim_config,
            )

            if enzyme_pdb is None or not enzyme_pdb.exists():
                logger.warning(
                    f"Cannot compute binding preference for '{cond.label}': "
                    f"enzyme PDB not found. Set enzyme_pdb_for_sasa in "
                    f"{settings_label} or contacts analysis settings."
                )
            else:
                bp_result = compute_condition_binding_preference(
                    cond=cond,
                    sim_config=sim_config,
                    analysis_dir=analysis_dir,
                    enzyme_pdb=enzyme_pdb,
                    threshold=settings["surface_exposure_threshold"],
                    include_default_aa_groups=settings["include_default_aa_groups"],
                    custom_protein_groups=settings["protein_groups"],
                    protein_partitions=settings["protein_partitions"],
                    polymer_type_selections=settings["polymer_type_selections"],
                )

            if bp_result is not None:
                logger.info(f"  Computed binding preference for {cond.label} at {temperature_K} K")
                base_data = {
                    "bp_result": bp_result,
                    "temperature_K": temperature_K,
                    "cond_label": cond.label,
                    "config_path": str(cond.config),
                }
                base_data.update(
                    self._binding_preference_extra_condition_data(cond, analysis_dir, bp_result)
                )
                return base_data

            logger.warning(
                f"Failed to compute binding preference for '{cond.label}'. "
                f"Ensure contacts_rep{{N}}.json files exist in {analysis_dir}."
            )
        else:
            logger.warning(
                f"No binding preference data found for '{cond.label}'. "
                f"Set compute_binding_preference=true or run contacts analysis "
                f"with binding preference enabled first."
            )

        base_data = {
            "bp_result": None,
            "temperature_K": temperature_K,
            "cond_label": cond.label,
            "config_path": str(cond.config),
        }
        base_data.update(
            self._binding_preference_extra_condition_data(cond, analysis_dir, bp_result)
        )
        return base_data

    def _find_contacts_analysis_dir(
        self,
        sim_config: Any,
        cond: "ConditionConfig",
        condition_output_dir: Path | None = None,
    ) -> Path:
        """Find the contacts analysis directory for a condition.

        Parameters
        ----------
        sim_config : Any
            Simulation configuration.
        cond : ConditionConfig
            Condition configuration.
        condition_output_dir : Path | None, optional
            Condition-specific contacts output directory.

        Returns
        -------
        Path
            Analysis directory path.
        """
        if condition_output_dir is not None:
            if condition_output_dir.exists():
                return condition_output_dir
            return condition_output_dir

        from polyzymd.compare.comparators._utils import find_analysis_dir

        return find_analysis_dir(
            sim_config,
            analysis_subdir="analysis/contacts",
            cond_config_path=Path(cond.config),
        )

    def _try_load_cached_binding_preference(
        self,
        cond: "ConditionConfig",
        analysis_dir: Path,
    ) -> Any:
        """Try loading cached binding preference for a condition.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        analysis_dir : Path
            Contacts analysis directory.

        Returns
        -------
        Any
            Cached binding preference result, or None if not found.
        """
        from polyzymd.compare.comparators._binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        return try_load_cached_binding_preference(cond, analysis_dir)

    def _resolve_compute_settings(self) -> dict[str, Any]:
        """Resolve compute settings with fallback to contacts settings.

        Returns
        -------
        dict[str, Any]
            Unified settings for binding preference computation.
        """
        settings = self.analysis_settings
        contacts_settings = None
        if hasattr(self.config, "analysis_settings"):
            contacts_settings = self.config.analysis_settings.get("contacts")

        def _get(attr: str, default: Any = None) -> Any:
            val = getattr(settings, attr, None)
            if val is not None:
                return val
            if contacts_settings is not None:
                val = getattr(contacts_settings, attr, None)
                if val is not None:
                    return val
            return default

        return {
            "enzyme_pdb_for_sasa": _get("enzyme_pdb_for_sasa"),
            "surface_exposure_threshold": _get("surface_exposure_threshold", 0.2),
            "include_default_aa_groups": _get("include_default_aa_groups", True),
            "protein_groups": _get("protein_groups"),
            "protein_partitions": _get("protein_partitions"),
            "polymer_type_selections": _get("polymer_type_selections"),
        }
