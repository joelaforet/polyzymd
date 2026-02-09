"""Test that all public modules can be imported."""

from pathlib import Path

import pytest


class TestImports:
    """Test basic package imports."""

    def test_import_polyzymd(self):
        """Test main package import."""
        import polyzymd

        assert hasattr(polyzymd, "__version__")

    def test_import_config(self):
        """Test config module import."""
        from polyzymd.config import SimulationConfig

        assert SimulationConfig is not None

    def test_import_schema(self):
        """Test schema imports."""
        from polyzymd.config.schema import (
            EnzymeConfig,
            OutputConfig,
            PolymerConfig,
            SimulationPhasesConfig,
            SolventConfig,
            SubstrateConfig,
            ThermodynamicsConfig,
        )

        assert EnzymeConfig is not None
        assert SubstrateConfig is not None
        assert PolymerConfig is not None
        assert SolventConfig is not None
        assert ThermodynamicsConfig is not None
        assert SimulationPhasesConfig is not None
        assert OutputConfig is not None

    def test_version_format(self):
        """Test version string format."""
        import polyzymd

        version = polyzymd.__version__
        # Should be semver format: X.Y.Z
        parts = version.split(".")
        assert len(parts) >= 2, f"Version {version} should have at least major.minor"
        assert parts[0].isdigit(), f"Major version should be numeric: {parts[0]}"
        assert parts[1].isdigit(), f"Minor version should be numeric: {parts[1]}"


class TestConfigValidation:
    """Test configuration validation."""

    def test_enzyme_config_required_fields(self):
        """Test EnzymeConfig requires name and pdb_path."""
        from pydantic import ValidationError

        from polyzymd.config.schema import EnzymeConfig

        # Should fail without required fields
        with pytest.raises(ValidationError):
            EnzymeConfig()

        # Should succeed with required fields
        config = EnzymeConfig(name="TestEnzyme", pdb_path="test.pdb")
        assert config.name == "TestEnzyme"
        assert config.pdb_path == Path("test.pdb")

    def test_monomer_probabilities_sum(self):
        """Test that monomer probabilities must sum to 1.0."""
        from pydantic import ValidationError

        from polyzymd.config.schema import MonomerSpec, PolymerConfig

        # Should fail if probabilities don't sum to 1.0
        with pytest.raises(ValidationError):
            PolymerConfig(
                type_prefix="TEST",
                monomers=[
                    MonomerSpec(label="A", probability=0.5, name="MonA"),
                    MonomerSpec(label="B", probability=0.3, name="MonB"),  # Sum = 0.8
                ],
                length=5,
                count=2,
                sdf_directory=Path("/tmp/test"),  # Required for cached mode
            )

        # Should succeed if probabilities sum to 1.0
        config = PolymerConfig(
            type_prefix="TEST",
            monomers=[
                MonomerSpec(label="A", probability=0.7, name="MonA"),
                MonomerSpec(label="B", probability=0.3, name="MonB"),
            ],
            length=5,
            count=2,
            sdf_directory=Path("/tmp/test"),  # Required for cached mode
        )
        assert len(config.monomers) == 2

    def test_thermodynamics_config(self):
        """Test ThermodynamicsConfig validation and defaults."""
        from pydantic import ValidationError

        from polyzymd.config.schema import ThermodynamicsConfig

        # Temperature is required (no default)
        with pytest.raises(ValidationError):
            ThermodynamicsConfig()

        # Should succeed with temperature provided
        config = ThermodynamicsConfig(temperature=300.0)
        assert config.temperature == 300.0
        assert config.pressure == 1.0  # Default pressure
