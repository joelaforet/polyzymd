"""An rc.5 config keeps working, unchanged, in every later 1.3 release.

PolyzyMD 1.3.0rc5 changed the config layout once; later 1.3 releases may only
add. ``schema_rc5.json`` is the contract of the rc.5 schema: for every key
reachable from ``SimulationConfig``, whether it is required, its default and
the values it accepts. A later schema must keep every key, require no key that
was optional, accept every value it accepted, and keep every default (a new
default would change the physics of an rc.5 config without an edit). The
example configs of rc.5 are kept as written and must still load.

Regenerate the contract only for a release that may break configs:
``python tests/config/test_schema_freeze.py --write``.
"""

from __future__ import annotations

import enum
import json
import sys
import typing
from pathlib import Path
from typing import Any

import pytest
from pydantic import BaseModel

from polyzymd.config.schema import SimulationConfig

DATA = Path(__file__).parent / "data"
CONTRACT = DATA / "schema_rc5.json"
CONFIGS = sorted((DATA / "rc5_configs").glob("*.yaml"))


def _values(annotation: Any) -> list[str] | None:
    """The values a Literal or Enum annotation accepts (inside Optional or a list)."""
    if isinstance(annotation, type) and issubclass(annotation, enum.Enum):
        return sorted(str(member.value) for member in annotation)
    if typing.get_origin(annotation) is typing.Literal:
        return sorted(str(value) for value in typing.get_args(annotation))
    for arg in typing.get_args(annotation):
        values = _values(arg)
        if values is not None:
            return values
    return None


def _models(annotation: Any) -> list[type[BaseModel]]:
    """The config models an annotation holds (a model, a list of them, or a union)."""
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        return [annotation]
    found: list[type[BaseModel]] = []
    for arg in typing.get_args(annotation):
        found += [m for m in _models(arg) if m not in found]
    return found


def _default(field: Any) -> Any:
    if field.is_required():
        return None
    value = field.get_default(call_default_factory=True)
    if isinstance(value, BaseModel):
        value = value.model_dump(mode="json")
    elif isinstance(value, enum.Enum):
        value = value.value
    return json.loads(json.dumps(value, default=str))


def contract(model: type[BaseModel] = SimulationConfig, path: str = "") -> dict[str, Any]:
    """Every key reachable from ``model``, by dotted path; a solute's keys under its type."""
    keys: dict[str, Any] = {}
    for name, field in model.model_fields.items():
        key = f"{path}{name}"
        keys[key] = {
            "required": field.is_required(),
            "default": _default(field),
            "values": _values(field.annotation),
        }
        for sub in _models(field.annotation):
            # A member of a union tagged by a one-value ``type`` (a solute) is keyed by its tag.
            tag = sub.model_fields.get("type")
            tags = _values(tag.annotation) if tag is not None else None
            branch = f"{key}[{tags[0]}]." if tags and len(tags) == 1 else f"{key}."
            keys.update(contract(sub, branch))
    return keys


def test_every_rc5_key_keeps_its_meaning() -> None:
    frozen = json.loads(CONTRACT.read_text())
    now = contract()
    problems = []
    for key, old in frozen.items():
        new = now.get(key)
        if new is None:
            problems.append(f"{key}: removed")
            continue
        if new["required"] and not old["required"]:
            problems.append(f"{key}: now required")
        if new["default"] != old["default"]:
            problems.append(f"{key}: default {old['default']!r} -> {new['default']!r}")
        if old["values"] is not None and not set(old["values"]) <= set(new["values"] or []):
            problems.append(
                f"{key}: no longer accepts {sorted(set(old['values']) - set(new['values'] or []))}"
            )
    added = [key for key in now if key not in frozen and now[key]["required"]]
    problems += [f"{key}: a new key that is required" for key in added]
    assert not problems, "an rc.5 config would need editing:\n" + "\n".join(problems)


@pytest.mark.parametrize("path", CONFIGS, ids=[p.stem for p in CONFIGS])
def test_the_rc5_example_configs_load_unchanged(path: Path) -> None:
    from polyzymd.config.loader import load_config

    config = load_config(path)
    assert config.solutes and config.solvent.water.model


def test_the_contract_covers_the_config() -> None:
    """A contract that missed whole sections would freeze nothing."""
    frozen = json.loads(CONTRACT.read_text())
    for key in (
        "solutes[protein].charges.method",
        "solutes[small_molecule].charges.model",
        "solutes[polymer].packing.tolerance",
        "solvent.water.model",
        "solvent.ions.nacl_concentration",
        "restraints.atom1.component",
        "modeling_exceptions.reason",
        "simulation_phases.production.time_step",
        "simulation_phases.minimization.freeze_solute",
    ):
        assert key in frozen, key


if __name__ == "__main__" and "--write" in sys.argv:
    CONTRACT.write_text(json.dumps(contract(), indent=1, sort_keys=True) + "\n")
    print(f"wrote {CONTRACT} ({len(contract())} keys)")
