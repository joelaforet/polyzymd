"""Regression test: all typed-result plugins must declare ReplicateResultClass.

This prevents a repeat of the HPC aggregate deserialization bug where
replicate results loaded from disk as raw dicts instead of Pydantic models.
"""

from __future__ import annotations

import pytest

from polyzymd.analyses import get_analysis, list_analyses

# All plugins that return typed Pydantic models from run_replicate
# and need ReplicateResultClass set for disk round-trip deserialization
TYPED_REPLICATE_PLUGINS = {
    "rmsf": "RMSFResult",
    "rmsd": "RMSDResult",
    "rg": "RgResult",
    "catalytic_triad": "TriadResult",
    "distances": "DistanceResult",
    "secondary_structure": "SecondaryStructureResult",
    "contacts": "ContactResult",
    "polymer_bridging": "PolymerBridgingReplicateResult",
    "sasa": "SASAResult",
    "hydrogen_bonds": "HydrogenBondResult",
}


class TestReplicateResultClassContract:
    """Every typed-result plugin must set ReplicateResultClass."""

    @pytest.mark.parametrize("plugin_name,expected_class_name", TYPED_REPLICATE_PLUGINS.items())
    def test_replicate_result_class_is_set(self, plugin_name, expected_class_name):
        """Plugin should have ReplicateResultClass pointing to the correct model."""
        cls = get_analysis(plugin_name)
        assert cls.ReplicateResultClass is not None, (
            f"{plugin_name} is missing ReplicateResultClass — "
            "HPC aggregate deserialization will return raw dicts instead of Pydantic models"
        )
        assert cls.ReplicateResultClass.__name__ == expected_class_name, (
            f"{plugin_name}.ReplicateResultClass should be {expected_class_name}, "
            f"got {cls.ReplicateResultClass.__name__}"
        )

    @pytest.mark.parametrize("plugin_name", TYPED_REPLICATE_PLUGINS.keys())
    def test_replicate_result_class_has_model_validate_json(self, plugin_name):
        """ReplicateResultClass must support Pydantic deserialization."""
        cls = get_analysis(plugin_name)
        result_cls = cls.ReplicateResultClass
        assert result_cls is not None
        assert hasattr(result_cls, "model_validate_json") or hasattr(
            result_cls, "load"
        ), f"{plugin_name}.ReplicateResultClass must have model_validate_json() or load()"

    def test_all_compute_stage_plugins_covered(self):
        """Sanity check for compute-stage plugin coverage."""
        all_plugins = list_analyses()
        for name, cls in all_plugins.items():
            if not cls.has_compute_stage:
                continue
            if name not in TYPED_REPLICATE_PLUGINS:
                # Plugin uses dict results so ReplicateResultClass should remain None
                assert cls.ReplicateResultClass is None, (
                    f"Plugin '{name}' has compute_stage and ReplicateResultClass set "
                    "but is not in TYPED_REPLICATE_PLUGINS — add it to the contract test"
                )
