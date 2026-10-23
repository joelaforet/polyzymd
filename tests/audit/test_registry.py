"""The check catalog: stable IDs, one implementation each, helpful lookups."""

from __future__ import annotations

import pytest

from polyzymd.audit import registry
from polyzymd.audit.registry import CheckSpec


@pytest.fixture
def empty_registry(monkeypatch):
    monkeypatch.setattr(registry, "_REGISTRY", {})
    monkeypatch.setattr(registry, "_IMPLEMENTATIONS", {})


def test_a_registered_check_is_found_with_its_implementation(empty_registry) -> None:
    spec = CheckSpec("geometry.example", "system", "FAIL", "an example mistake")

    @registry.register(spec)
    def run(*args):
        return "ran"

    assert registry.get("geometry.example") is spec
    assert registry.implementation("geometry.example")() == "ran"
    assert registry.catalog() == [spec]


def test_an_id_cannot_be_registered_twice(empty_registry) -> None:
    spec = CheckSpec("parameters.example", "component", "FAIL", "x")
    registry.register(spec)(lambda: None)
    with pytest.raises(ValueError, match="registered twice"):
        registry.register(spec)


def test_an_unknown_id_names_the_closest(empty_registry) -> None:
    registry.register(CheckSpec("geometry.buried_water", "system", "WARN", "x"))(lambda: None)
    with pytest.raises(KeyError, match="did you mean geometry.buried_water"):
        registry.get("geometry.buried_waters")


def test_a_check_can_be_limited_to_component_types() -> None:
    spec = CheckSpec("composition.polymer", "component", "INFO", "x", component_types=("polymer",))
    assert spec.applies_to("polymer")
    assert not spec.applies_to("protein")
