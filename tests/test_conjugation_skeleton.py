"""Tests for Phase 0-1 conjugation skeleton."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation import CovalentModificationBuilder
from polyzymd.builders.conjugation.exceptions import ConjugationNotImplementedError
from polyzymd.config.schema import (
    ConjugationChainPolicyConfig,
    ConjugationConfig,
    SimulationConfig,
)


def _minimal_simulation_config_data() -> dict:
    """Build a minimal valid simulation config dictionary.

    Returns
    -------
    dict
        Dictionary suitable for :meth:`SimulationConfig.model_validate`.
    """
    return {
        "name": "conjugation-skeleton",
        "enzyme": {"name": "enzyme", "pdb_path": "enzyme.pdb"},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {
                    "name": "eq",
                    "duration": 0.1,
                    "samples": 10,
                    "ensemble": "NVT",
                    "temperature": 300.0,
                }
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
                "time_step": 2.0,
            },
        },
    }


class TestConjugationConfigParsing:
    """Configuration parsing tests for the conjugation skeleton."""

    def test_simulation_config_accepts_absent_conjugation(self):
        """Top-level conjugation remains optional for legacy configs."""
        config = SimulationConfig.model_validate(_minimal_simulation_config_data())
        assert config.conjugation is None

    def test_simulation_config_accepts_disabled_conjugation(self):
        """Explicit disabled conjugation parses without changing defaults."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {"enabled": False}

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.enabled is False
        assert config.conjugation.chain_policy.protein_chain == "A"
        assert config.conjugation.chain_policy.moiety_chain == "C"

    def test_enabled_ingest_existing_block_parses(self):
        """Enabled ingest_existing config parses with Pablo policy placeholders."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "mode": "ingest_existing",
            "source_pdb_path": "prebuilt_conjugate.pdb",
            "ccd_pablo": {"lookup_policy": "offline_cached"},
        }

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.enabled is True
        assert config.conjugation.mode.value == "ingest_existing"
        assert config.conjugation.source_pdb_path.name == "prebuilt_conjugate.pdb"

    def test_enabled_construct_block_parses(self):
        """Enabled construct config parses attachment placeholders."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "mode": "construct",
            "attachments": [
                {
                    "name": "lys23-peg",
                    "site": {"chain_id": "A", "residue_name": "LYS", "residue_number": 23},
                    "moiety": {"name": "PEG", "smiles": "COCCO"},
                    "mechanism": {"name": "amide"},
                }
            ],
        }

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.mode.value == "construct"
        assert config.conjugation.attachments[0].moiety.name == "PEG"

    def test_invalid_mode_rejected(self):
        """Unknown conjugation modes fail validation."""
        with pytest.raises(ValidationError):
            ConjugationConfig(enabled=True, mode="unsupported")

    def test_invalid_ccd_policy_rejected(self):
        """Unknown CCD/Pablo policies fail validation."""
        with pytest.raises(ValidationError):
            ConjugationConfig(
                enabled=True,
                ccd_pablo={"lookup_policy": "always_guess"},
            )

    def test_chain_policy_defaults(self):
        """Chain policy defaults preserve the PolyzyMD convention."""
        policy = ConjugationChainPolicyConfig()
        assert policy.protein_chain == "A"
        assert policy.substrate_chain == "B"
        assert policy.moiety_chain == "C"
        assert policy.solvent_start_chain == "D"

    def test_chain_policy_rejects_overlap(self):
        """Component and solvent chain IDs must not overlap."""
        with pytest.raises(ValidationError, match="must not overlap"):
            ConjugationChainPolicyConfig(solvent_start_chain="C")

    def test_chain_policy_rejects_nonalpha(self):
        """Chain IDs must be alphabetic single characters."""
        with pytest.raises(ValidationError, match="alphabetic"):
            ConjugationChainPolicyConfig(protein_chain="1")


class TestCovalentModificationBuilder:
    """Builder behavior tests for the no-op skeleton."""

    def test_from_config_accepts_simulation_config(self):
        """Factory extracts the nested conjugation config."""
        sim_config = MagicMock()
        sim_config.conjugation = ConjugationConfig(enabled=False)

        builder = CovalentModificationBuilder.from_config(sim_config)

        assert builder.config is sim_config.conjugation

    def test_disabled_builder_returns_input_unchanged(self):
        """Disabled conjugation is a no-op with empty sidecar models."""
        topology = object()
        builder = CovalentModificationBuilder(ConjugationConfig(enabled=False))

        result = builder.build(topology)

        assert result.topology is topology
        assert result.metadata.enabled is False
        assert result.metadata.components == []
        assert result.diagnostics.enabled is False
        assert result.diagnostics.diagnostics == []

    def test_enabled_construct_raises_clear_not_implemented(self):
        """Construct mode fails clearly in the Phase 0-1 skeleton."""
        builder = CovalentModificationBuilder(ConjugationConfig(enabled=True, mode="construct"))

        with pytest.raises(ConjugationNotImplementedError, match="not implemented"):
            builder.build(object())

    def test_enabled_ingest_existing_raises_clear_not_implemented(self):
        """Ingest mode routes through the Pablo adapter boundary and fails clearly."""
        builder = CovalentModificationBuilder(
            ConjugationConfig(
                enabled=True,
                mode="ingest_existing",
                source_pdb_path="prebuilt_conjugate.pdb",
            )
        )

        with pytest.raises(ConjugationNotImplementedError, match="Pablo|not implemented"):
            builder.build(object())


class TestSystemBuilderConjugationHook:
    """SystemBuilder compatibility tests for the pre-solvation hook."""

    def test_absent_conjugation_returns_combined_topology(self):
        """Legacy configs without conjugation remain unchanged."""
        from polyzymd.builders.system_builder import SystemBuilder

        topology = object()
        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = topology

        config = MagicMock()
        config.conjugation = None

        assert builder._apply_conjugation(config) is topology

    def test_disabled_conjugation_returns_combined_topology(self):
        """Disabled conjugation remains a no-op through SystemBuilder."""
        from polyzymd.builders.system_builder import SystemBuilder

        topology = object()
        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = topology

        config = MagicMock()
        config.conjugation = ConjugationConfig(enabled=False)

        assert builder._apply_conjugation(config) is topology

    def test_enabled_conjugation_fails_before_solvation(self):
        """Enabled skeleton paths fail clearly at the pre-solvation hook."""
        from polyzymd.builders.system_builder import SystemBuilder

        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = object()
        builder._enzyme_topology = object()
        builder._working_dir = None
        builder._n_enzyme_molecules = 1
        builder._n_substrate_molecules = 0
        builder._n_polymer_chains = 0

        config = MagicMock()
        config.conjugation = ConjugationConfig(enabled=True, mode="construct")

        with pytest.raises(ConjugationNotImplementedError, match="not implemented"):
            builder._apply_conjugation(config)
