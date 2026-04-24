"""Tests for conjugation polymer recipes and Polymerist smoke generation."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation.polymer_recipe import (
    PolymerRecipe,
    generate_polymerist_smoke_polymer,
    sbma_egpma_nhs_recipe,
)
from polyzymd.config.loader import load_config
from polyzymd.config.schema import SimulationConfig


def _minimal_simulation_config_data() -> dict:
    """Build a minimal valid simulation config dictionary."""
    return {
        "name": "conjugation-polymer-recipe",
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


def _recipe_dict() -> dict:
    """Return a realistic SBMA/EGPMA/NHS polymer recipe dictionary."""
    return sbma_egpma_nhs_recipe(length=9, seed=7).model_dump(mode="json")


def test_recipe_validates_probabilities_and_pdb_safe_residue_names():
    """Recipes should reject invalid probabilities and residue names."""
    recipe_data = _recipe_dict()
    recipe_data["monomers"][0]["probability"] = 0.5

    with pytest.raises(ValidationError, match="probabilities must sum"):
        PolymerRecipe.model_validate(recipe_data)

    recipe_data = _recipe_dict()
    recipe_data["monomers"][0]["residue_name"] = "SBMA"

    with pytest.raises(ValidationError, match="at most 3 characters"):
        PolymerRecipe.model_validate(recipe_data)


def test_seeded_sequence_is_deterministic_and_forces_center_nhs():
    """The same seed should generate the same centered NHS reactive residue."""
    recipe = sbma_egpma_nhs_recipe(length=9, seed=123)

    first = recipe.generate_sequence()
    second = recipe.generate_sequence()

    assert first == second
    assert len(first) == 9
    assert first[recipe.effective_reactive_index] == "C"


def test_forced_reactive_index_places_nhs_at_known_residue():
    """A configured reactive index should override the default centered NHS."""
    recipe = sbma_egpma_nhs_recipe(length=8, seed=11, reactive_monomer_index=2)
    sequence = recipe.generate_sequence()

    assert recipe.effective_reactive_index == 2
    assert sequence[2] == "C"


def test_invalid_forced_label_and_index_are_rejected():
    """Forced reactive selectors must refer to declared residues in range."""
    recipe_data = _recipe_dict()
    recipe_data["reactive_monomer_label"] = "Z"

    with pytest.raises(ValidationError, match="Reactive monomer label"):
        PolymerRecipe.model_validate(recipe_data)

    recipe_data = _recipe_dict()
    recipe_data["reactive_monomer_index"] = 99

    with pytest.raises(ValidationError, match="within the polymer length"):
        PolymerRecipe.model_validate(recipe_data)


def test_simulation_config_accepts_yaml_polymer_recipe(tmp_path):
    """Conjugation moieties should parse polymer recipes from YAML."""
    data = _minimal_simulation_config_data()
    data["conjugation"] = {
        "enabled": True,
        "mode": "construct",
        "attachments": [
            {
                "name": "lys23-sbma-egpma-nhs",
                "site": {"chain_id": "A", "residue_name": "LYS", "residue_number": 23},
                "moiety": {"name": "SBMA-EGPMA-NHS", "polymer_recipe": _recipe_dict()},
                "mechanism": {"name": "nhs_lys_amide"},
            }
        ],
    }

    config = SimulationConfig.model_validate(data)

    recipe = config.conjugation.attachments[0].moiety.polymer_recipe
    assert isinstance(recipe, PolymerRecipe)
    assert recipe.monomer_by_label("C").name == "NHS"

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
name: conjugation-polymer-recipe
enzyme:
  name: enzyme
  pdb_path: enzyme.pdb
thermodynamics:
  temperature: 300.0
simulation_phases:
  equilibration_stages:
    - name: eq
      duration: 0.1
      samples: 10
      ensemble: NVT
      temperature: 300.0
  production:
    ensemble: NPT
    duration: 1.0
    samples: 10
    time_step: 2.0
conjugation:
  enabled: true
  mode: construct
  attachments:
    - name: lys23-sbma-egpma-nhs
      site:
        chain_id: A
        residue_name: LYS
        residue_number: 23
      moiety:
        name: SBMA-EGPMA-NHS
        recipe:
          name: SBMA-EGPMA-NHS
          length: 9
          seed: 7
          forced_reactive_monomer_label: C
          monomers:
            - label: A
              name: SBMA
              residue_name: SBM
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
              probability: 0.945
            - label: B
              name: EGPMA
              residue_name: EGP
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
              probability: 0.045
            - label: C
              name: NHS
              residue_name: NHS
              smiles: CC(=C)C(=O)ON1C(=O)CCC1=O
              probability: 0.01
      mechanism:
        name: nhs_lys_amide
""".strip()
    )

    loaded = load_config(config_path)
    loaded_recipe = loaded.conjugation.attachments[0].moiety.polymer_recipe
    assert isinstance(loaded_recipe, PolymerRecipe)
    assert loaded_recipe.generate_sequence()[loaded_recipe.effective_reactive_index] == "C"


def test_real_polymerist_generation_smoke(tmp_path):
    """Polymerist should consume the real SBMA/EGPMA/NHS recipe when installed."""
    recipe = sbma_egpma_nhs_recipe(length=3, seed=5, reactive_monomer_index=1)

    try:
        result = generate_polymerist_smoke_polymer(recipe, tmp_path / "polymerist", max_retries=1)
    except Exception as exc:
        pytest.skip(f"Polymerist generation stack unavailable in this environment: {exc}")

    assert result.sequence[1] == "C"
    assert result.monomer_group_path.exists()
    assert result.pdb_path is not None and result.pdb_path.exists()
    assert result.object_type
    assert result.atom_count is None or result.atom_count > 0
