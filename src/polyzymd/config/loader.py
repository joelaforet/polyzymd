"""
YAML configuration loader and saver for PolyzyMD.

This module provides functions to load and save SimulationConfig
objects from/to YAML files, with support for Path objects and
environment variable expansion.
"""

from __future__ import annotations

import difflib
import logging
import os
import re
from collections.abc import Hashable
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Union, get_args

import click
import yaml
from pydantic import BaseModel, ValidationError

from polyzymd.config.schema import SimulationConfig
from polyzymd.core.branding import prepend_file_header

#: Config keys whose values are file or directory paths, resolved against the
#: folder holding the config when they are relative.
PATH_KEYS = frozenset(
    {
        "pdb_path",
        "custom_substructures_path",
        "sdf_path",
        "sdf_directory",
        "cache_directory",
        "initiation",
        "polymerization",
        "termination",
    }
)
#: Where runs and job files go. They are resolved against the config's folder,
#: as the input files are, so a command finds the runs from any folder; they
#: are never copied or hashed as inputs.
OUTPUT_PATH_KEYS = frozenset({"projects_directory", "scratch_directory"})
LOGGER = logging.getLogger(__name__)

#: A ``$VAR`` or ``${VAR}`` left in a path after expansion: the variable is unset.
UNSET_VARIABLE = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)")


class ConfigFileError(click.ClickException, ValueError):
    """A config file that cannot be read: not YAML, not a mapping, or a key given twice.

    It is a ``ValueError`` for the code that catches those, and a Click error,
    so every command prints it as ``error:`` and ``fix:`` lines and exits 1.
    """

    def __init__(self, path: Path, problem: str, fix: str) -> None:
        super().__init__(f"{path}: {problem}")
        self.fix = fix

    def __str__(self) -> str:
        return f"{self.message}\nfix: {self.fix}"

    def show(self, file: Any = None) -> None:
        """Print the ``error:`` and ``fix:`` lines."""
        click.echo(f"error: {self.message}", err=True)
        click.echo(f"fix: {self.fix}", err=True)


class _UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that refuses a key given twice in one mapping."""

    def construct_mapping(self, node: yaml.MappingNode, deep: bool = False) -> Dict[Any, Any]:
        # Keys merged in with "<<:" may be overridden, so only the mapping's own
        # keys, which flatten_mapping puts after the merged ones, are checked.
        own = sum(1 for key_node, _ in node.value if key_node.tag != "tag:yaml.org,2002:merge")
        self.flatten_mapping(node)
        seen: set[Any] = set()
        for key_node, _ in node.value[len(node.value) - own :]:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, Hashable):
                break  # the parent reports the unhashable key
            if key in seen:
                raise yaml.constructor.ConstructorError(
                    None, None, f"duplicate key {key!r}", key_node.start_mark
                )
            seen.add(key)
        return super().construct_mapping(node, deep=deep)


def read_yaml_mapping(path: Path) -> Dict[str, Any]:
    """Read a YAML file that must hold a mapping, refusing duplicate keys.

    Raises
    ------
    ConfigFileError
        When the file is not valid YAML, holds a key twice, or is not a mapping.
    """
    try:
        with open(path, "r") as f:
            data = yaml.load(f, Loader=_UniqueKeyLoader)
    except yaml.YAMLError as error:
        mark = getattr(error, "problem_mark", None)
        where = f" at line {mark.line + 1}, column {mark.column + 1}" if mark else ""
        problem = getattr(error, "problem", None) or str(error)
        raise ConfigFileError(
            path,
            f"not valid YAML: {problem}{where}",
            "Correct the YAML at that line (indent with spaces, not tabs; give each key once).",
        ) from error
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigFileError(
            path,
            f"must be a mapping of config sections, not a {type(data).__name__}",
            "Write the config as 'key: value' sections (name:, engine:, enzyme:, ...).",
        )
    return data


def _expand_paths(data: Dict[str, Any], base_path: Path) -> Dict[str, Any]:
    """Recursively expand relative paths in configuration data.

    Converts relative paths to absolute paths based on the config file location.
    Also expands ``~`` and environment variables in path strings. A path that
    names an unset variable is kept as written; new builds refuse it
    (:meth:`SimulationConfig.require_buildable`).

    Sentinel values (e.g. ``"default"``) are passed through untouched so that
    downstream Pydantic validators can resolve them to bundled resources.

    Args:
        data: Configuration dictionary
        base_path: Directory containing the config file

    Returns:
        Configuration with expanded paths
    """
    path_keys = PATH_KEYS | OUTPUT_PATH_KEYS

    # Sentinel values that should be forwarded to Pydantic validators as-is,
    # not treated as filesystem paths.
    _SENTINEL_VALUES = {"default"}

    def expand_value(key: str, value: Any) -> Any:
        if key in path_keys and isinstance(value, str):
            # Pass through sentinel values without path expansion
            if value.lower().strip() in _SENTINEL_VALUES:
                return value
            expanded = os.path.expanduser(os.path.expandvars(value))
            if UNSET_VARIABLE.search(expanded):
                return value
            path = Path(expanded)
            # Convert relative paths to absolute based on config file location
            if not path.is_absolute():
                path = base_path / path
            # Earlier versions did not expand "~" and wrote the runs under a
            # folder named "~" beside the config; keep finding those runs.
            old = base_path / value
            if key in OUTPUT_PATH_KEYS and value.startswith("~") and not path.exists():
                if old.is_dir():
                    LOGGER.warning(
                        "%s %r: using %s, where earlier versions put the runs; '~' now means "
                        "the home folder. Move the runs and the config finds them there.",
                        key,
                        value,
                        old,
                    )
                    return str(old)
            return str(path)
        elif isinstance(value, dict):
            return {k: expand_value(k, v) for k, v in value.items()}
        elif isinstance(value, list):
            return [expand_value(key, item) for item in value]
        return value

    return {k: expand_value(k, v) for k, v in data.items()}


def _convert_paths_to_relative(data: Dict[str, Any], base_path: Path) -> Dict[str, Any]:
    """Convert absolute paths to relative paths for saving.

    Args:
        data: Configuration dictionary with absolute paths
        base_path: Directory where config file will be saved

    Returns:
        Configuration with relative paths
    """
    path_keys = PATH_KEYS | OUTPUT_PATH_KEYS

    # Sentinel values that should be forwarded as-is (see _expand_paths).
    _SENTINEL_VALUES = {"default"}

    def relativize_value(key: str, value: Any) -> Any:
        if key in path_keys and isinstance(value, str):
            # Pass through sentinel values without relativizing
            if value.lower().strip() in _SENTINEL_VALUES:
                return value
            path = Path(value)
            if path.is_absolute():
                try:
                    return str(path.relative_to(base_path))
                except ValueError:
                    # Path is not relative to base_path, keep absolute
                    return value
            return value
        elif isinstance(value, dict):
            return {k: relativize_value(k, v) for k, v in value.items()}
        elif isinstance(value, list):
            return [relativize_value(key, item) for item in value]
        return value

    return {k: relativize_value(k, v) for k, v in data.items()}


class ConfigLoader:
    """Custom YAML loader with support for includes and references."""

    def __init__(self, base_path: Path):
        self.base_path = base_path

    def load(self, path: Path) -> Dict[str, Any]:
        """Read the YAML mapping at ``path`` and expand its paths."""
        return _expand_paths(read_yaml_mapping(path), self.base_path)


#: The fix for a config that lacks a whole required section.
TEMPLATE_FIX = (
    "start from the complete commented config that polyzymd study add-condition LABEL --new "
    "writes in a study folder (polyzymd study init DIR makes one), and copy your values into it."
)
#: How YAML names the Python types a wrong value can have.
_YAML_TYPE_NAMES = {
    str: "text",
    int: "a number",
    float: "a number",
    bool: "true/false",
    list: "a list",
    dict: "a mapping",
    type(None): "empty",
}


def _models_in(annotation: Any) -> list[type[BaseModel]]:
    """Return the config models an annotation holds (``X``, ``X | None``, ``list[X]`` ...)."""
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        return [annotation]
    return [model for arg in get_args(annotation) for model in _models_in(arg)]


@lru_cache(maxsize=1)
def _schema_paths() -> tuple[tuple[str, ...], ...]:
    """Return the path of every field of :class:`SimulationConfig`, nested ones included."""
    paths: list[tuple[str, ...]] = []

    def walk(model: type[BaseModel], prefix: tuple[str, ...], seen: tuple[type, ...]) -> None:
        for name, field in model.model_fields.items():
            paths.append((*prefix, name))
            for sub in _models_in(field.annotation):
                if sub not in seen:
                    walk(sub, (*prefix, name), (*seen, sub))

    walk(SimulationConfig, (), (SimulationConfig,))
    return tuple(paths)


def _annotation_at(keys: tuple[str, ...]) -> Any:
    """Return the annotation of the field at ``keys``, or None when it is not a field."""
    annotation: Any = SimulationConfig
    for key in keys:
        for model in _models_in(annotation):
            if key in model.model_fields:
                annotation = model.model_fields[key].annotation
                break
        else:
            return None
    return annotation


def _shape(annotation: Any) -> str | None:
    """Describe the YAML a field expects: a mapping or a list of mappings, with its keys."""
    models = _models_in(annotation)
    if len(models) != 1:
        return None
    names = list(models[0].model_fields)
    keys = " and ".join([", ".join(names[:-1]), names[-1]]) if len(names) > 1 else names[0]
    many = getattr(annotation, "__origin__", None) is list
    return f"a list of mappings with {keys}" if many else f"a mapping with {keys}"


def _key_score(key: str, name: str) -> float:
    """Score how likely ``name`` is the key a user meant by ``key`` (0 to 1)."""
    key, name = key.lower().replace("-", "_"), name.lower()
    score = difflib.SequenceMatcher(None, key, name).ratio()
    # A key that is a whole word of a short field name: pdb -> pdb_path.
    short, long = sorted((key, name), key=len)
    if short in long.split("_") and 3 * len(short) >= len(long):
        score = max(score, 0.8)
    return score


def _closest_path(parent: tuple[str, ...], key: str) -> tuple[str, ...] | None:
    """Return the field path closest to the unknown ``key`` under ``parent``, or None."""
    best, best_rank = None, None
    for path in _schema_paths():
        score = _key_score(key, path[-1])
        if score < 0.75:
            continue
        rank = (score + (0.05 if path[:-1] == parent else 0.0), -len(path))
        if best_rank is None or rank > best_rank:
            best, best_rank = path, rank
    return best


def _unknown_key_messages(loc: tuple[Any, ...], value: Any) -> list[str]:
    """Explain an unknown key, with the closest field path when there is one.

    A mapping under an unknown key (``system: {pdb: ...}``) is explained key
    by key, so each value is pointed to the field it belongs in.
    """
    where = ".".join(str(part) for part in loc)
    parent = tuple(part for part in loc[:-1] if isinstance(part, str))

    def mean(path: tuple[str, ...]) -> str:
        # A field at the same level keeps the list indices of the location.
        if path[:-1] == parent:
            return ".".join(str(part) for part in (*loc[:-1], path[-1]))
        return ".".join(path)

    if isinstance(value, dict) and value:
        inner = {str(k): _closest_path((), str(k)) for k in value}
        if any(inner.values()):
            return [
                f"{where}.{k}: unknown key" + (f"; did you mean {mean(p)}?" if p else "")
                for k, p in inner.items()
            ]
    path = _closest_path(parent, str(loc[-1]))
    if path is not None:
        return [f"{where}: unknown key; did you mean {mean(path)}?"]
    valid = [p[-1] for p in _schema_paths() if p[:-1] == parent]
    return [
        (
            f"{where}: unknown key; valid keys here: {', '.join(valid)}"
            if valid
            else f"{where}: unknown key"
        )
    ]


def config_error_lines(error: ValidationError, path: Union[str, Path]) -> list[tuple[str, str]]:
    """Explain a config ``ValidationError`` as ``error:`` and ``fix:`` lines.

    Unknown keys name the closest valid field path, at that level or in
    another section; a missing section points to the commented template;
    a value of the wrong type names the shape the field expects. Other
    errors keep the message of the schema.

    Parameters
    ----------
    error : ValidationError
        The error raised while validating the config.
    path : str or Path
        The config file, for the fix line.

    Returns
    -------
    list of tuple of str
        ``("error", text)`` and ``("fix", text)`` pairs, in printing order.
    """
    lines: list[tuple[str, str]] = []
    fixes: list[str] = []
    ours = error.title == SimulationConfig.__name__
    for item in error.errors():
        loc = tuple(item["loc"])
        where = ".".join(str(part) for part in loc) or "config"
        keys = tuple(part for part in loc if isinstance(part, str))
        kind = item["type"]
        annotation = _annotation_at(keys) if ours and keys else None
        shape = _shape(annotation) if annotation is not None else None
        if ours and kind == "extra_forbidden" and loc:
            lines += [("error", text) for text in _unknown_key_messages(loc, item.get("input"))]
        elif ours and kind == "missing":
            what = "required section missing" if shape else "required key missing"
            lines.append(("error", f"{where}: {what}" + (f"; expected {shape}" if shape else "")))
            if len(keys) == 1 and shape and TEMPLATE_FIX not in fixes:
                fixes.append(TEMPLATE_FIX)
        elif shape and kind in ("model_type", "dict_type", "model_attributes_type", "list_type"):
            got = _YAML_TYPE_NAMES.get(type(item.get("input")), type(item.get("input")).__name__)
            lines.append(("error", f"{where}: expected {shape}, not {got}"))
        else:
            message, _, fix = item["msg"].partition("\nfix: ")
            lines.append(("error", f"{where}: {message}"))
            if fix and fix not in fixes:
                fixes.append(fix)
    fixes.append(f"correct these keys in {path} and run polyzymd validate again.")
    return lines + [("fix", text) for text in fixes]


def load_config(path: Union[str, Path]) -> SimulationConfig:
    """Load a SimulationConfig from a YAML file.

    Parameters
    ----------
    path : str or Path
        Path to the YAML configuration file.

    Returns
    -------
    SimulationConfig
        Validated configuration instance.

    Raises
    ------
    FileNotFoundError
        If the config file doesn't exist.
    ConfigFileError
        If the file is not YAML, is not a mapping, or gives a key twice.
    pydantic.ValidationError
        If the configuration is invalid.

    Examples
    --------
    >>> config = load_config("my_simulation.yaml")
    >>> print(config.enzyme.name)
    "LipA"
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    base_path = path.parent.absolute()

    data = ConfigLoader(base_path).load(path)

    return SimulationConfig.model_validate(data)


def save_config(
    config: SimulationConfig, path: Union[str, Path], relative_paths: bool = True
) -> None:
    """Save a SimulationConfig to a YAML file.

    Parameters
    ----------
    config : SimulationConfig
        Configuration to save.
    path : str or Path
        Destination path for the YAML file.
    relative_paths : bool, optional
        Whether to convert paths to relative, by default True.

    Examples
    --------
    >>> config = SimulationConfig(...)
    >>> save_config(config, "output_config.yaml")
    """
    path = Path(path)

    # Create parent directory if needed
    path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to dict, handling Path objects
    data = config.model_dump(mode="json")

    if relative_paths:
        data = _convert_paths_to_relative(data, path.parent.absolute())

    header = prepend_file_header("", comment_prefix="#")

    # Custom YAML representer for cleaner output — use a local Dumper
    # subclass so we don't mutate the global yaml.Dumper state.
    class _CleanDumper(yaml.Dumper):
        pass

    def str_representer(dumper: yaml.Dumper, data: str) -> yaml.Node:
        if "\n" in data:
            return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    _CleanDumper.add_representer(str, str_representer)

    with open(path, "w") as f:
        f.write(header)
        yaml.dump(
            data,
            f,
            Dumper=_CleanDumper,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
            width=100,
        )
