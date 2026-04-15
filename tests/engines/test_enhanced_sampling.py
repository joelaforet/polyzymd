"""Tests for enhanced-sampling extension hooks on SimulationEngine.

Verifies the contract established by Phase 7:

- Optional hooks are NOT abstract — engines can be instantiated without them.
- Default implementations raise ``NotImplementedError`` with engine-specific
  messages.
- Both ``OpenMMEngine`` and ``GromacsEngine`` inherit the correct messages.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine
from polyzymd.engines.bias import EnhancedSamplingProtocol, ExternalBiasSpec

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture()
def bias_spec() -> ExternalBiasSpec:
    """Minimal PLUMED bias specification."""
    return ExternalBiasSpec(
        bias_type="plumed",
        script_content="d: DISTANCE ATOMS=1,2\nPRINT ARG=d FILE=dist.dat",
    )


@pytest.fixture()
def protocol() -> EnhancedSamplingProtocol:
    """Minimal metadynamics protocol definition."""
    return EnhancedSamplingProtocol(
        protocol_type="metadynamics",
        parameters={"hill_height": 1.0, "hill_width": 0.2, "pace": 500},
    )


@pytest.fixture()
def submit_request(tmp_path: Path) -> EngineSubmitRequest:
    """Minimal submission request for hook tests."""
    return EngineSubmitRequest(
        replicate=1,
        config_path=tmp_path / "config.yaml",
        working_dir=tmp_path / "run_1",
    )


# ============================================================================
# ABC contract tests
# ============================================================================


class TestHooksAreNotAbstract:
    """Enhanced-sampling hooks must be optional — not abstract."""

    def test_attach_external_bias_is_not_abstract(self) -> None:
        """attach_external_bias should NOT be in __abstractmethods__."""
        assert "attach_external_bias" not in SimulationEngine.__abstractmethods__

    def test_configure_enhanced_sampling_is_not_abstract(self) -> None:
        """configure_enhanced_sampling should NOT be in __abstractmethods__."""
        assert "configure_enhanced_sampling" not in SimulationEngine.__abstractmethods__

    def test_save_engine_state_is_not_abstract(self) -> None:
        """save_engine_state should NOT be in __abstractmethods__."""
        assert "save_engine_state" not in SimulationEngine.__abstractmethods__

    def test_save_bias_state_is_not_abstract(self) -> None:
        """save_bias_state should NOT be in __abstractmethods__."""
        assert "save_bias_state" not in SimulationEngine.__abstractmethods__


# ============================================================================
# Default error message tests — OpenMM
# ============================================================================


class TestOpenMMEnhancedSamplingHooks:
    """OpenMMEngine should raise engine-specific NotImplementedError."""

    @pytest.fixture()
    def engine(self) -> object:
        """Create an OpenMMEngine with minimal mock config."""
        from polyzymd.engines.openmm.engine import OpenMMEngine

        config = MagicMock()
        return OpenMMEngine(config=config)

    def test_attach_external_bias_raises(
        self, engine: object, submit_request: EngineSubmitRequest, bias_spec: ExternalBiasSpec
    ) -> None:
        """OpenMM should raise with engine name in message."""
        with pytest.raises(NotImplementedError, match="not yet supported for openmm"):
            engine.attach_external_bias(submit_request, bias_spec)

    def test_configure_enhanced_sampling_raises(
        self, engine: object, protocol: EnhancedSamplingProtocol
    ) -> None:
        """OpenMM should raise with engine name and future-release note."""
        with pytest.raises(NotImplementedError, match="not yet supported for openmm"):
            engine.configure_enhanced_sampling(protocol)

    def test_configure_enhanced_sampling_mentions_future_release(
        self, engine: object, protocol: EnhancedSamplingProtocol
    ) -> None:
        """Error message should guide contributors to the extension point."""
        with pytest.raises(NotImplementedError, match="future release"):
            engine.configure_enhanced_sampling(protocol)

    def test_save_engine_state_raises(self, engine: object, tmp_path: Path) -> None:
        """save_engine_state should raise for OpenMM."""
        with pytest.raises(NotImplementedError, match="not yet supported for openmm"):
            engine.save_engine_state(tmp_path, replicate=1)

    def test_save_bias_state_raises(self, engine: object, tmp_path: Path) -> None:
        """save_bias_state should raise for OpenMM."""
        with pytest.raises(NotImplementedError, match="not yet supported for openmm"):
            engine.save_bias_state(tmp_path, replicate=1)


# ============================================================================
# Default error message tests — GROMACS
# ============================================================================


class TestGromacsEnhancedSamplingHooks:
    """GromacsEngine should raise engine-specific NotImplementedError."""

    @pytest.fixture()
    def engine(self) -> object:
        """Create a GromacsEngine with minimal mock config."""
        from polyzymd.engines.gromacs.engine import GromacsEngine

        config = MagicMock()
        config.gromacs = MagicMock()
        config.gromacs.gmx_binary = "gmx"
        return GromacsEngine(config=config, gmx_binary="gmx")

    def test_attach_external_bias_raises(
        self, engine: object, submit_request: EngineSubmitRequest, bias_spec: ExternalBiasSpec
    ) -> None:
        """GROMACS should raise with engine name in message."""
        with pytest.raises(NotImplementedError, match="not yet supported for gromacs"):
            engine.attach_external_bias(submit_request, bias_spec)

    def test_configure_enhanced_sampling_raises(
        self, engine: object, protocol: EnhancedSamplingProtocol
    ) -> None:
        """GROMACS should raise with engine name and future-release note."""
        with pytest.raises(NotImplementedError, match="not yet supported for gromacs"):
            engine.configure_enhanced_sampling(protocol)

    def test_configure_enhanced_sampling_mentions_future_release(
        self, engine: object, protocol: EnhancedSamplingProtocol
    ) -> None:
        """Error message should guide contributors to the extension point."""
        with pytest.raises(NotImplementedError, match="future release"):
            engine.configure_enhanced_sampling(protocol)

    def test_save_engine_state_raises(self, engine: object, tmp_path: Path) -> None:
        """save_engine_state should raise for GROMACS."""
        with pytest.raises(NotImplementedError, match="not yet supported for gromacs"):
            engine.save_engine_state(tmp_path, replicate=1)

    def test_save_bias_state_raises(self, engine: object, tmp_path: Path) -> None:
        """save_bias_state should raise for GROMACS."""
        with pytest.raises(NotImplementedError, match="not yet supported for gromacs"):
            engine.save_bias_state(tmp_path, replicate=1)


# ============================================================================
# Bias model tests
# ============================================================================


class TestBiasModels:
    """Verify ExternalBiasSpec and EnhancedSamplingProtocol Pydantic models."""

    def test_external_bias_spec_defaults(self) -> None:
        """ExternalBiasSpec should default to plumed with no script."""
        spec = ExternalBiasSpec()
        assert spec.bias_type == "plumed"
        assert spec.script_path is None
        assert spec.script_content is None

    def test_external_bias_spec_with_script_path(self, tmp_path: Path) -> None:
        """ExternalBiasSpec should accept a file path."""
        plumed_file = tmp_path / "plumed.dat"
        plumed_file.write_text("PRINT ARG=d FILE=dist.dat")
        spec = ExternalBiasSpec(bias_type="plumed", script_path=plumed_file)
        assert spec.script_path == plumed_file

    def test_external_bias_spec_with_inline_script(self) -> None:
        """ExternalBiasSpec should accept inline script content."""
        spec = ExternalBiasSpec(script_content="d: DISTANCE ATOMS=1,2")
        assert "DISTANCE" in spec.script_content

    def test_enhanced_sampling_protocol_creation(self) -> None:
        """EnhancedSamplingProtocol should store protocol type and params."""
        proto = EnhancedSamplingProtocol(
            protocol_type="metadynamics",
            parameters={"hill_height": 1.0, "pace": 500},
        )
        assert proto.protocol_type == "metadynamics"
        assert proto.parameters["hill_height"] == 1.0

    def test_enhanced_sampling_protocol_empty_params(self) -> None:
        """EnhancedSamplingProtocol should allow empty parameters."""
        proto = EnhancedSamplingProtocol(protocol_type="umbrella")
        assert proto.parameters == {}
