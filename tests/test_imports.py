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


class TestCoSolventVolumeValidation:
    """Test co-solvent volume fraction / concentration validation."""

    def test_concentration_cosolvent_does_not_crash_validator(self):
        """Concentration-based co-solvents must not crash validate_volume_fractions.

        Regression: sum() over volume_fraction raised TypeError when any
        co-solvent used concentration (volume_fraction=None).
        """
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[CoSolventSpec(name="urea", concentration=2.0)]
        )
        assert len(config.co_solvents) == 1

    def test_mixed_volume_fraction_and_concentration_cosolvents(self):
        """Mixed volume_fraction + concentration co-solvents should validate."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[
                CoSolventSpec(name="dmso", volume_fraction=0.3),
                CoSolventSpec(name="urea", concentration=2.0),
            ]
        )
        assert len(config.co_solvents) == 2

    def test_volume_fractions_exceeding_one_rejected(self):
        """Total volume fractions >= 1.0 should still be rejected."""
        from pydantic import ValidationError

        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        with pytest.raises(ValidationError, match="volume fraction must be < 1.0"):
            SolventConfig(
                co_solvents=[
                    CoSolventSpec(name="dmso", volume_fraction=0.6),
                    CoSolventSpec(name="ethanol", volume_fraction=0.5),
                ]
            )

    def test_volume_fractions_below_one_accepted(self):
        """Total volume fractions < 1.0 should pass validation."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[
                CoSolventSpec(name="dmso", volume_fraction=0.3),
                CoSolventSpec(name="ethanol", volume_fraction=0.2),
            ]
        )
        assert len(config.co_solvents) == 2


# ---------------------------------------------------------------------------
# B10 – save_config must not mutate global yaml.Dumper
# ---------------------------------------------------------------------------


class TestSaveConfigYamlDumper:
    """save_config uses a local Dumper subclass, not the global yaml.Dumper."""

    def test_global_dumper_unmodified_after_save(self, tmp_path):
        """Calling save_config must not add representers to yaml.Dumper."""
        import yaml

        from polyzymd.config.loader import save_config
        from polyzymd.config.schema import SimulationConfig

        # Snapshot the global Dumper's str representer before save_config
        original_representer = yaml.Dumper.yaml_representers.get(str)

        # Build a minimal config — save_config needs a valid SimulationConfig.
        # Use model_construct to skip validation (we only care about YAML dump).
        config = SimulationConfig.model_construct(
            name="test",
            _fields_set=set(),
        )
        out = tmp_path / "test_config.yaml"
        save_config(config, out)

        after_representer = yaml.Dumper.yaml_representers.get(str)
        assert after_representer is original_representer, (
            "save_config must not mutate yaml.Dumper.yaml_representers"
        )

    def test_multiline_strings_use_block_style(self, tmp_path):
        """Multiline strings in saved YAML use literal block (|) style."""
        import yaml

        from polyzymd.config.loader import save_config
        from polyzymd.config.schema import SimulationConfig

        config = SimulationConfig.model_construct(
            name="multi\nline\ntest",
            _fields_set=set(),
        )
        out = tmp_path / "test_config.yaml"
        save_config(config, out)

        content = out.read_text()
        # The literal block indicator (| or |-) should appear for multiline strings
        assert "|-" in content or "|\n" in content


# ---------------------------------------------------------------------------
# B12 – reaction template paths included in path resolution
# ---------------------------------------------------------------------------


class TestReactionPathResolution:
    """_expand_paths and _convert_paths_to_relative handle reaction template keys."""

    def test_expand_paths_resolves_initiation(self):
        """Relative initiation path should be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"initiation": "templates/init.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["initiation"] == "/configs/templates/init.rxn"

    def test_expand_paths_resolves_polymerization(self):
        """Relative polymerization path should be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"polymerization": "templates/poly.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["polymerization"] == "/configs/templates/poly.rxn"

    def test_expand_paths_resolves_termination(self):
        """Relative termination path should be expanded to absolute."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"termination": "templates/term.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["termination"] == "/configs/templates/term.rxn"

    def test_expand_paths_preserves_absolute_reaction_path(self):
        """Absolute reaction paths should not be modified."""
        from polyzymd.config.loader import _expand_paths

        data = {"reactions": {"initiation": "/absolute/init.rxn"}}
        base = Path("/configs")
        result = _expand_paths(data, base)
        assert result["reactions"]["initiation"] == "/absolute/init.rxn"

    def test_convert_paths_relativizes_reaction_paths(self):
        """Absolute reaction paths should be converted to relative for saving."""
        from polyzymd.config.loader import _convert_paths_to_relative

        data = {"reactions": {"initiation": "/configs/templates/init.rxn"}}
        base = Path("/configs")
        result = _convert_paths_to_relative(data, base)
        assert result["reactions"]["initiation"] == "templates/init.rxn"
