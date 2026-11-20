"""Create a study folder: the layout ``polyzymd study init`` writes.

See the "Study folders" explanation page for what each part is for. Each
condition lives in ``conditions/<name>/``: either a template config to fill
in, or a copy of an existing ``config.yaml`` with the input files it names
copied into ``conditions/<name>/structures/`` and its paths made relative,
so the folder can be moved and published as a whole. Either way, the
config writes its runs into the git-ignored ``runs/`` folder of the
project, or of the study when it is in no project (:func:`runs_folder`).
"""

from __future__ import annotations

import os
import re
import shutil
import textwrap
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Mapping, Sequence

from polyzymd.analyses.exceptions import ProtocolError

#: Folders of a new study, with what each holds.
FOLDERS = {
    "conditions": "One folder per condition: its config.yaml and structures/",
    "structures": "Reference structures for analyses, such as a crystal structure",
    "analyses": "Your measurement functions, listed in study.yaml as file.py:function",
    "figures": "Notebooks and scripts that turn stored results into the paper's figures",
    "results": "Stored per-replicate values, reports and figures, written by polyzymd analyze",
    "environment": "The software environment that produced the results",
}

GITIGNORE = """\
# Where this machine keeps the trajectories: never commit or publish it.
data.local.yaml
# The runs of the conditions' configs, unless they set scratch_directory.
runs/
# Fragments and chains that dynamic polymer builds write into the folder they run in.
.polymer_cache/
# Python and notebook caches.
__pycache__/
*.pyc
.ipynb_checkpoints/
# What polyzymd study freeze lays out for upload.
deposit/
# Full logs of polyzymd commands; the console shows only warnings.
logs/
# Job scripts and logs of polyzymd analyze --submit, which name this
# machine's paths, and SLURM logs of the simulations.
results/*/slurm/
conditions/*/slurm_logs/
# Software environments, which are rebuilt, never committed.
.pixi/
.venv/
"""

DATA_EXAMPLE = """\
# data.local.yaml says where THIS machine keeps each condition's run
# directories. It is gitignored and never published, because a path is only
# a pointer to where the data lives now: moving data must not change the study.
#
# Without data.local.yaml, each condition's config.yaml says where its runs are
# (its scratch_directory), which suits whoever ran the simulations. After
# downloading trajectories, for example from Zenodo, write it with
#
#     polyzymd study locate DOWNLOAD_DIR
#
# or by hand, as condition label -> directory holding its run directories:
#
# No polymer: /path/to/downloaded/no_polymer
# SBMA 50%: /pl/active/my_lab/polyzymd_sims/LipA_363K
"""

MIT = """\
MIT License

Copyright (c) {year} {holder}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

CC_BY = """\
Copyright (c) {year} {holder}

The data, results and figures of this study (everything except the code in
analyses/ and figures/, which LICENSE-code covers) are licensed under the
Creative Commons Attribution 4.0 International License (CC BY 4.0):
https://creativecommons.org/licenses/by/4.0/legalcode

SPDX-License-Identifier: CC-BY-4.0

Replace this file and LICENSE-code to publish under other licences, and set
metadata.license in study.yaml to match.
"""


@dataclass
class CreatedStudy:
    """What :func:`create_study` wrote."""

    root: Path
    conditions: dict[str, Path]
    copied: dict[str, list[str]] = field(default_factory=dict)
    left_absolute: dict[str, list[str]] = field(default_factory=dict)
    commit: str | None = None


def condition_folder(label: str) -> str:
    """Return the folder name of a condition label: lower case, words joined by ``_``."""
    return re.sub(r"[^\w.+-]+", "_", label.strip().lower()).strip("_") or "condition"


#: Folder names a study or condition may not take: polyzymd writes folders of
#: these names, and leaves job, log and run folders out of git and the deposit.
RESERVED_FOLDERS = (
    *("runs", "logs", "slurm", "slurm_logs", "deposit", "results"),
    *("analyses", "stats", "figures", "conditions", "structures", "environment"),
)


def check_label(label: str, what: str) -> None:
    """Refuse an empty study or condition ``label``, or one whose folder name is reserved.

    Reserved folder names are those in :data:`RESERVED_FOLDERS`. ``what``
    (``"study"`` or ``"condition"``) names the label in the message.
    """
    if not label.strip():
        raise ProtocolError(f"The {what} label is empty.", hint=f"Give the {what} a name.")
    name = condition_folder(label)
    if name in RESERVED_FOLDERS:
        raise ProtocolError(
            f"The {what} label {label!r} gives the folder name {name!r}, which is reserved: "
            "polyzymd uses folders of that name.",
            hint=f"Choose another label. Reserved: {', '.join(RESERVED_FOLDERS)}.",
        )


#: What add-condition prints, and copied configs say, about where the runs go.
RUNS_WARNING = (
    "the runs go into {runs} unless you set scratch_directory in config.yaml; trajectories "
    "can use a lot of disk space, so on a cluster set scratch_directory to scratch storage"
)


#: The header of a copied config whose runs exist: the copy keeps them.
REUSE_NOTE = (
    "scratch_directory is the folder of the runs of {config}: this config reuses them, "
    "and its new replicates go there too"
)

#: The header of a copied config that keeps the absolute scratch_directory of the original.
KEPT_NOTE = (
    "scratch_directory is kept from {config}: the runs go there, into the same "
    "run folders as those of {config}"
)

#: The header of a new condition copied with --from whose original sets an absolute scratch_directory.
OWN_NOTE = (
    "scratch_directory is a folder of this condition's own, in the scratch storage of "
    "the {config} it was copied from, so the two never share run folders"
)


def runs_folder(study_root: Path, label: str) -> Path:
    """Return the folder that a new condition's config writes its runs into.

    It is ``runs/<study>/<condition>`` in the project folder when the folder
    above the study holds a ``project.yaml``, and ``runs/<condition>`` in
    the study folder otherwise. Both ``runs/`` folders are git-ignored.
    """
    root = Path(study_root).resolve()
    home = runs_home(root)
    return home / "runs" / root.relative_to(home) / condition_folder(label)


def runs_home(study_root: Path) -> Path:
    """Return the folder whose ``runs/`` holds the study's runs: its project's, or its own."""
    from polyzymd.analyses.project_file import PROJECT_FILE

    root = Path(study_root).resolve()
    return root.parent if (root.parent / PROJECT_FILE).is_file() else root


def ignore_runs(folder: Path) -> None:
    """Add ``runs/`` to ``folder/.gitignore`` when the file exists without it.

    Studies and projects made before the runs went into ``runs/`` lack it.
    """
    gitignore = Path(folder) / ".gitignore"
    if gitignore.is_file():
        lines = gitignore.read_text().splitlines()
        if "runs/" not in lines:
            gitignore.write_text("\n".join([*lines, "runs/", ""]))


def copy_condition(
    config: Path,
    folder: Path,
    runs: Path,
    *,
    name: str | None = None,
    own_folder: str | None = None,
    source_folder: str | None = None,
) -> tuple[list[str], list[str]]:
    """Copy ``config`` to ``folder/config.yaml`` with the input files it names.

    Every path key of the config (``config.loader.PATH_KEYS``) that names an
    existing file or directory is copied into ``folder/structures/`` and
    written as a path relative to the new config. ``projects_directory``
    becomes ``runs`` (the folder :func:`runs_folder` gives), relative to the
    new config. ``scratch_directory`` becomes the folder that holds the
    config's runs, when it has any (:func:`existing_runs`), so the copy
    reuses them: analysis, ``status``, ``run`` and ``submit`` of either
    config see the same replicates, and new ones go beside them. Otherwise
    an absolute ``scratch_directory`` (or one starting with an environment
    variable) is kept as written, and a relative or missing one becomes
    null, so new runs go into ``runs``. With ``own_folder``, the copy is to
    become another condition (``add-condition --from``) and gets runs of its
    own: existing runs are not reused, and a kept ``scratch_directory`` gets
    ``own_folder`` appended, after removing ``source_folder`` from its end
    when the original was itself such a copy. The run folders of the two
    configs, named by the same ``naming_template``, then never coincide.
    A polymer ``cache_directory`` is
    left out so the default is used. The config hash leaves these out, so
    stored results still match. The text is kept with its comments and
    indentation, and only the lines of the copied paths, the output
    directories and the name change. A config whose changed lines cannot be
    found (a path written in flow style) is written again, with list items
    indented under their key and without comments.
    With ``name``, the top-level ``name:`` of the copy is set to it, as
    :func:`new_condition` sets the template's name from the folder.

    Returns
    -------
    tuple of list of str
        The copied paths and the paths left absolute because nothing exists
        there.
    """
    import yaml

    from polyzymd.config.loader import PATH_KEYS, _expand_paths
    from polyzymd.config.schema import SimulationConfig

    config = Path(config).expanduser().resolve()
    try:
        raw = yaml.safe_load(config.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(f"Cannot read {config}: {exc}", hint="Give a config.yaml.") from exc
    data = _expand_paths(raw, config.parent)
    structures = folder / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    missing: list[str] = []
    seen: dict[Path, str] = {}

    def place(source: Path) -> str:
        if source in seen:
            return seen[source]
        name, stem, n = source.name, source.stem, 1
        while (structures / name).exists():
            n += 1
            name = f"{stem}_{n}{source.suffix}"
        target = structures / name
        if source.is_dir():
            shutil.copytree(source, target)
        else:
            shutil.copy2(source, target)
        seen[source] = f"structures/{name}"
        copied.append(str(source))
        return seen[source]

    def walk(key: str | None, value: Any) -> Any:
        if isinstance(value, dict):
            return {k: walk(k, v) for k, v in value.items()}
        if isinstance(value, list):
            return [walk(key, item) for item in value]
        if key in PATH_KEYS and isinstance(value, str) and value.lower().strip() != "default":
            source = Path(value)
            if key != "cache_directory" and source.exists():
                return place(source)
            if key != "cache_directory":
                missing.append(value)
        return value

    rewritten = walk(None, data)
    polymers = rewritten.get("polymers")
    dropped_cache = isinstance(polymers, dict) and "cache_directory" in polymers
    if isinstance(polymers, dict):
        polymers.pop("cache_directory", None)
    existing = existing_runs(config) if own_folder is None else None
    scratch = _absolute_scratch(raw.get("output")) if existing is None else None
    if scratch is not None and own_folder is not None:
        scratch = scratch.rstrip("/")
        if source_folder and scratch.endswith(f"/{source_folder}"):
            scratch = scratch[: -len(source_folder) - 1]
        scratch = f"{scratch}/{own_folder}"
    output = rewritten.setdefault("output", {})
    if isinstance(output, dict):
        output["projects_directory"] = os.path.relpath(runs, folder.resolve())
        output["scratch_directory"] = _from(folder, existing) if existing is not None else scratch
    if name is not None:
        rewritten["name"] = name
    target = folder / "config.yaml"
    # Keep the original text, comments and indentation included: change only
    # the lines of the copied paths, the output folders and the name.
    lines = config.read_text().splitlines(keepends=True)
    while lines and lines[0].startswith(("# Copied by polyzymd", "#   ")):
        lines.pop(0)  # the header of an earlier copy
    kept = _with_path_lines("".join(lines), _changed_paths(raw, rewritten, PATH_KEYS))
    if dropped_cache:
        kept = re.sub(r"^[ \t]+cache_directory:[^\n]*\n", "", kept, count=1, flags=re.MULTILINE)
    kept = _with_output_lines(kept, rewritten.get("output"))
    if name is not None:
        kept = _with_name_line(kept, name)
    try:
        same = yaml.safe_load(kept) == rewritten
    except yaml.YAMLError:
        same = False
    if not same:
        kept = yaml.dump(rewritten, Dumper=_indented_list_dumper(), sort_keys=False)
    if existing is not None:
        note = REUSE_NOTE.format(config=config.name)
    elif scratch is not None:
        note = (KEPT_NOTE if own_folder is None else OWN_NOTE).format(config=config.name)
    else:
        where = f"{os.path.relpath(runs, folder.resolve())} (relative to this file)"
        note = RUNS_WARNING.format(runs=where)
    target.write_text(
        f"# Copied by polyzymd from {config.name}\n"
        + "".join(f"#   {line}\n" for line in textwrap.wrap(f"{note}.", 74))
        + kept
    )
    try:
        SimulationConfig.from_yaml(target)
    except (OSError, ValueError) as exc:
        raise ProtocolError(
            f"The copy of {config} in {target} does not load: {exc}",
            hint="Check that the original config loads with polyzymd validate.",
        ) from exc
    return copied, missing


#: The header line of a config that ``--set`` had to write again, without its comments.
COMMENTS_LOST_NOTE = "--set could not keep the comments of this config"


def read_assignments(entries: Sequence[str], option: str) -> list[tuple[str, Any]]:
    """Parse ``KEY=VALUE`` entries of ``option``, reading each value as YAML.

    Raises
    ------
    ProtocolError
        If an entry has no ``=`` or no key, or its value is not YAML.
    """
    import yaml

    pairs: list[tuple[str, Any]] = []
    for entry in entries:
        key, separator, value = entry.partition("=")
        key = key.strip()
        if not separator or not key:
            raise ProtocolError(
                f"Cannot read {option} {entry!r}.",
                hint=f"Write it as {option} KEY=VALUE, for example "
                f"{option} thermodynamics.temperature=330.",
            )
        try:
            pairs.append((key, yaml.safe_load(value)))
        except yaml.YAMLError as exc:
            raise ProtocolError(
                f"Cannot read the value of {option} {key!r}: {exc}",
                hint=f"Quote the value, for example {option} description='SBMA 50%'.",
            ) from exc
    return pairs


def key_parts(key: str) -> list[str | int]:
    """Split a dotted config key such as ``solutes[1].monomers[0].probability`` into parts."""
    parts: list[str | int] = []
    for segment in key.split("."):
        match = re.fullmatch(r"([^\[\]\s]+)((?:\[\d+\])*)", segment.strip())
        if match is None:
            raise ProtocolError(
                f"Cannot read the config key {key!r}.",
                hint="Write a dotted path with list indices, for example "
                "solutes[1].monomers[0].probability.",
            )
        parts.append(match.group(1))
        parts += [int(index) for index in re.findall(r"\d+", match.group(2))]
    return parts


def _key_text(parts: Sequence[str | int]) -> str:
    return "".join(f"[{p}]" if isinstance(p, int) else f".{p}" for p in parts).lstrip(".")


def _set_value(data: dict, parts: list[str | int], value: Any) -> None:
    """Set ``parts`` of ``data`` to ``value``, adding missing mappings; never add list items."""
    from polyzymd.config.loader import _closest_path

    node: Any = data
    for i, part in enumerate(parts):
        where, above = _key_text(parts[: i + 1]), _key_text(parts[:i]) or "the config"
        if isinstance(part, int) and isinstance(node, list) and part >= len(node):
            raise ProtocolError(
                f"{where} is out of range: {above} has {len(node)} item(s), counted from 0.",
                hint="Use the index of an existing item; --set does not add list items.",
            )
        if not isinstance(node, list if isinstance(part, int) else dict):
            keys = tuple(p for p in parts[:i] if isinstance(p, str))
            path = _closest_path(keys[:-1], keys[-1]) if keys else None
            raise ProtocolError(
                f"{above} is not a {'list' if isinstance(part, int) else 'mapping'} in the "
                f"config, so {where} cannot be set.",
                hint=f"Did you mean {'.'.join(path)}?" if path else "Check the key.",
            )
        if i == len(parts) - 1:
            node[part] = value
        else:
            if isinstance(part, str) and node.get(part) is None and isinstance(parts[i + 1], str):
                node[part] = {}
            node = node.get(part) if isinstance(node, dict) else node[part]


def flow_yaml(value: Any) -> str:
    """Return ``value`` as one line of flow-style YAML."""
    import yaml

    text = yaml.safe_dump(value, default_flow_style=True, width=float("inf"), allow_unicode=True)
    return text.removesuffix("\n...\n").strip()


def _set_in_text(text: str, parts: list[str | int], value: Any) -> str | None:
    """Return ``text`` with ``parts`` set to ``value`` by one text edit, or ``None``.

    A scalar is replaced where it is written; a missing key goes into its
    block mapping, as ``key: value`` or ``key: {rest: value}``. Anything
    else (a whole list or section, a flow mapping) returns ``None``.
    """
    import yaml

    node = yaml.compose(text)
    rest = list(parts)
    while rest:
        part = rest[0]
        if isinstance(part, int) and isinstance(node, yaml.SequenceNode):
            node = node.value[part]
        elif isinstance(part, str) and isinstance(node, yaml.MappingNode):
            found = [v for k, v in node.value if getattr(k, "value", None) == part]
            if not found:
                break
            node = found[0]
        else:
            return None
        rest.pop(0)
    if not rest:
        if not isinstance(node, yaml.ScalarNode):
            return None
        return text[: node.start_mark.index] + flow_yaml(value) + text[node.end_mark.index :]
    if not isinstance(node, yaml.MappingNode) or node.flow_style or not node.value:
        return None
    if any(isinstance(part, int) for part in rest):
        return None
    for part in reversed(rest[1:]):
        value = {part: value}
    key, first = node.value[0]
    start = text.rfind("\n", 0, key.start_mark.index) + 1
    indent = text[start : key.start_mark.index]
    line = f"{' ' * key.start_mark.column}{flow_yaml(rest[0])}: {flow_yaml(value)}\n"
    if not indent.strip():
        return text[:start] + line + text[start:]
    if isinstance(first, yaml.ScalarNode) and first.end_mark.line == key.start_mark.line:
        end = text.find("\n", first.end_mark.index)
        end = len(text) if end < 0 else end + 1
        return text[:end] + ("" if text[:end].endswith("\n") else "\n") + line + text[end:]
    return None


def set_config_values(config: Path, settings: Sequence[tuple[str, Any]]) -> None:
    """Set each dotted key of ``settings`` in ``config``, keeping its comments.

    Each setting is one text edit (:func:`_set_in_text`), checked by reading
    the file back. When an edit cannot be made that way, the config is
    written again without its comments, keeping its header and adding
    :data:`COMMENTS_LOST_NOTE` to it.

    Raises
    ------
    ProtocolError
        If a key cannot be read, indexes past the end of a list, or goes
        through a value that is not a mapping or list.
    """
    import yaml

    text = config.read_text()
    data = yaml.safe_load(text) or {}
    kept = True
    for key, value in settings:
        parts = key_parts(key)
        _set_value(data, parts, value)
        if kept:
            edited = _set_in_text(text, parts, value)
            try:
                kept = edited is not None and yaml.safe_load(edited) == data
            except yaml.YAMLError:
                kept = False
            text = edited if kept else text
    if not kept:
        header = re.match(r"(?:#[^\n]*\n)*", text).group(0)
        text = (
            header
            + f"#   {COMMENTS_LOST_NOTE}.\n"
            + yaml.dump(data, Dumper=_indented_list_dumper(), sort_keys=False)
        )
    config.write_text(text)


def existing_runs(config: Path) -> Path | None:
    """Return the folder holding the runs of ``config``, or ``None`` when it has none.

    A relative folder is relative to the config's folder.
    """
    from polyzymd.config.schema import SimulationConfig

    try:
        loaded = SimulationConfig.from_yaml(Path(config))
        if not loaded.discover_replicate_dirs():
            return None
        where = Path(loaded.output.effective_scratch_directory).expanduser()
    except (OSError, ValueError):
        return None
    return (where if where.is_absolute() else Path(config).resolve().parent / where).resolve()


def _absolute_scratch(output: Any) -> str | None:
    """Return the ``scratch_directory`` of a raw ``output`` block when it is absolute.

    One starting with ``~`` or an environment variable counts as absolute:
    it names a machine's storage, not a folder beside the config.
    """
    value = output.get("scratch_directory") if isinstance(output, dict) else None
    if not isinstance(value, str) or not value.strip():
        return None
    value = value.strip()
    expanded = os.path.expanduser(os.path.expandvars(value))
    return value if value.startswith("$") or Path(expanded).is_absolute() else None


def _from(folder: Path, target: Path) -> str:
    """Return ``target`` relative to ``folder`` when they share a folder below the top level."""
    folder = folder.resolve()
    if len(Path(os.path.commonpath([folder, target])).parts) > 2:
        return os.path.relpath(target, folder)
    return str(target)


def _changed_paths(old: Any, new: Any, keys: frozenset[str], key: str | None = None) -> list:
    """List ``(key, old, new)`` for each path value of ``keys`` that differs, in file order."""
    if isinstance(old, dict) and isinstance(new, dict):
        return [c for k, v in old.items() if k in new for c in _changed_paths(v, new[k], keys, k)]
    if isinstance(old, list) and isinstance(new, list):
        return [c for a, b in zip(old, new) for c in _changed_paths(a, b, keys, key)]
    if key in keys and isinstance(old, str) and isinstance(new, str) and old != new:
        return [(key, old, new)]
    return []


def _with_path_lines(text: str, changes: list) -> str:
    """Return ``text`` with each ``key: old`` line of ``changes`` set to the new path.

    A value written in flow style (``{pdb_path: x}``) is not found; the caller
    then writes the whole config instead.
    """
    import json

    for key, old, new in changes:
        quoted = "|".join(re.escape(form) for form in (json.dumps(old), f"'{old}'", old))
        text = re.sub(
            rf"^([ \t]*(?:-[ \t]+)?{re.escape(key)}:[ \t]*)(?:{quoted})([ \t]+#.*)?[ \t]*$",
            lambda match, new=new: f"{match.group(1)}{json.dumps(new)}{match.group(2) or ''}",
            text,
            count=1,
            flags=re.MULTILINE,
        )
    return text


def _indented_list_dumper() -> type:
    import yaml

    class Dumper(yaml.SafeDumper):
        """Indent list items under their key, as the examples and tutorials do."""

        def increase_indent(self, flow: bool = False, indentless: bool = False) -> None:
            return super().increase_indent(flow, False)

    return Dumper


def _with_name_line(text: str, name: str) -> str:
    """Return ``text`` with its top-level ``name:`` set to ``name``, keeping the line's comment."""
    import json

    return re.sub(
        r"^name:[^#\n]*?(\s+#.*)?$",
        lambda match: f"name: {json.dumps(name)}{match.group(1) or ''}",
        text,
        count=1,
        flags=re.MULTILINE,
    )


def _with_output_lines(text: str, output: Any) -> str:
    """Return ``text`` with its ``projects_directory`` and ``scratch_directory`` set from ``output``.

    Only those lines change; their comments stay, except the note of a
    deposited config that its machine path was removed.
    """
    import json

    if not isinstance(output, dict):
        return text

    def value(match: re.Match) -> str:
        new = output.get(match.group(2))
        written = "null" if new is None else json.dumps(str(new))
        comment = match.group(3) or ""
        if "machine path removed by polyzymd" in comment:
            comment = ""
        return f"{match.group(1)}{match.group(2)}: {written}{comment}"

    text = re.sub(
        r"^(\s+)(projects_directory|scratch_directory):[^#\n]*?(\s+#.*)?$",
        value,
        text,
        flags=re.MULTILINE,
    )
    # A key the config leaves out (scratch_directory, usually) goes under output:.
    missing = [
        key
        for key in ("projects_directory", "scratch_directory")
        if key in output and not re.search(rf"^\s+{key}:", text, flags=re.MULTILINE)
    ]
    if not missing:
        return text
    block = re.search(r"^output:[ \t]*(#.*)?\n(?:([ \t]+)\S)?", text, flags=re.MULTILINE)
    if block is None:
        text = text if text.endswith("\n") or not text else text + "\n"
        text += "output:\n"
        at, indent = len(text), "  "
    else:
        at, indent = text.index("\n", block.start(0)) + 1, block.group(2) or "  "
    added = "".join(
        f"{indent}{key}: "
        + ("null" if output[key] is None else json.dumps(str(output[key])))
        + "\n"
        for key in missing
    )
    return text[:at] + added + text[at:]


def record_data_location(study_root: Path, label: str, config: Path) -> None:
    """Write where the runs of ``config`` are into the study's ``data.local.yaml``, for ``label``.

    The directory is the one holding the config's runs (:func:`existing_runs`).
    Nothing is written when the config has no run: the copy's runs then go
    into ``runs/``. Entries for other conditions are kept.
    """
    import yaml

    from polyzymd.analyses.study_file import DATA_FILE

    where = existing_runs(config)
    if where is None:
        return
    file = Path(study_root) / DATA_FILE
    current = yaml.safe_load(file.read_text()) if file.is_file() else None
    current = dict(current or {})
    current[label] = str(where.resolve())
    file.write_text(
        "# Where this machine keeps each condition's runs; never committed or published.\n"
        + yaml.safe_dump(current, sort_keys=False)
    )


def _study_yaml(conditions: dict[str, str], equilibration: str | None) -> str:
    import yaml

    listed = (
        yaml.safe_dump({"conditions": conditions}, sort_keys=False)
        if conditions
        else ("# Add each condition: label: conditions/<name>/config.yaml\nconditions: {}\n")
    )
    window = equilibration or "0ns"
    note = "" if equilibration else "   # set the window to discard from every replicate"
    return f"""\
# The analysis protocol of this study: see
# https://polyzymd.readthedocs.io/en/latest/how_to/study_yaml.html
equilibration: {window}{note}
{listed}analyses: {{}}
# analyses:
#   rg: {{}}
#   contacts:
#     method: occlusion
#   lid_opening:
#     function: analyses/lid.py:lid_distance
#     kind: timeseries
#     unit: A
#     selections:
#       lid: "protein and resid 140-150 and name CA"
#       core: "protein and resid 4-120 and name CA"
"""


def _readme(name: str) -> str:
    from polyzymd.citation import citation_line

    rows = "\n".join(f"| `{folder}/` | {what} |" for folder, what in FOLDERS.items())
    return f"""\
# {name}

A molecular dynamics study made with PolyzyMD. `study.yaml` holds the
analysis protocol: the conditions, the equilibration window and every
analysis setting. Each condition's simulation is in `conditions/<name>/config.yaml`.

| Path | Holds |
|---|---|
| `study.yaml` | The analysis protocol |
{rows}
| `data.example.yaml` | How to say where this machine keeps the trajectories (`data.local.yaml`) |
| `LICENSE-data`, `LICENSE-code` | The licences of the data and of the code |

## Reproduce

Install the PolyzyMD version the results were made with (each report and
`manifest.json` record it; see `environment/`), then:

1. **Figures, from the stored results only:** run the notebooks and scripts in
   `figures/`, which read `pz.Study("study.yaml").results(run)`.
2. **Analyses, from the trajectories:** download them, run
   `polyzymd study locate DOWNLOAD_DIR`, check with `polyzymd study check`,
   and run `polyzymd analyze --study study.yaml`.
3. **Simulations:** each `conditions/<name>/config.yaml` builds and runs its
   condition with `polyzymd`; the replicate number is the random seed.
   Results agree within the statistical noise of MD, not bit for bit.

## Publish

Fill in `metadata:` in `study.yaml` and run `polyzymd study freeze`, which
writes `manifest.json`, `CITATION.cff`, `.zenodo.json` and
`md_checklist.yaml`, commits and tags the study, and lays out `deposit/` for
upload.

## How to cite

Cite this study's paper, and PolyzyMD, which produced the analyses:

> {citation_line()}
"""


def _environment_readme() -> str:
    import polyzymd

    return f"""\
# Environment

This study was started with PolyzyMD {polyzymd.__version__}. To reproduce
it, install the version its results were made with, for example from
https://github.com/joelaforet/polyzymd at that tag, with its pixi environment:

    pixi install -e analysis

Every report records the PolyzyMD version that made it, and
`polyzymd study freeze` records PolyzyMD and its main dependencies in
`manifest.json` and warns when stored results were made with another version.
To pin the whole environment, copy the `pixi.toml` and `pixi.lock` you ran
with into this folder and commit them.
"""


def create_study(
    root: str | Path,
    *,
    conditions: dict[str, Path] | None = None,
    new_conditions: list[str] | None = None,
    equilibration: str | None = None,
    holder: str | None = None,
    git: bool = True,
) -> CreatedStudy:
    """Write a new study folder at ``root``.

    ``conditions`` maps labels to existing configs, copied with their input
    files (:func:`copy_condition`); ``new_conditions`` are labels for which a
    template config is written (:func:`new_condition`). With ``git``, the folder is made a
    git repository and everything in it is committed; nothing is committed
    afterwards.

    Raises
    ------
    ProtocolError
        If ``root`` already holds a ``study.yaml``, a label is given twice or
        is reserved (:func:`check_label`), or a config cannot be copied.
    """
    from polyzymd.analyses.study_file import STUDY_FILE
    from polyzymd.analyses.study_git import init_repository

    root = Path(root).expanduser().resolve()
    if (root / STUDY_FILE).exists():
        raise ProtocolError(
            f"{root} already holds a {STUDY_FILE}.",
            hint="Choose a new folder, or edit the existing study.yaml.",
        )
    labels = [*(conditions or {}), *(new_conditions or [])]
    for label in labels:
        check_label(label, "condition")
    folders = [condition_folder(label) for label in labels]
    if len(set(labels)) != len(labels) or len(set(folders)) != len(folders):
        raise ProtocolError(
            "Each condition needs its own label and folder name.",
            hint=f"The labels {labels} give the folders {folders}; rename the clashing ones.",
        )
    root.mkdir(parents=True, exist_ok=True)
    for folder in FOLDERS:
        (root / folder).mkdir(exist_ok=True)
    created = CreatedStudy(root, {})
    for label, config in (conditions or {}).items():
        folder = root / "conditions" / condition_folder(label)
        created.copied[label], created.left_absolute[label] = copy_condition(
            Path(config), folder, runs_folder(root, label)
        )
        record_data_location(root, label, Path(config))
        created.conditions[label] = folder / "config.yaml"
    for label in new_conditions or []:
        folder = root / "conditions" / condition_folder(label)
        new_condition(folder, runs_folder(root, label))
        created.conditions[label] = folder / "config.yaml"

    year, holder = date.today().year, holder or "the study's authors"
    listed = {label: str(path.relative_to(root)) for label, path in created.conditions.items()}
    (root / STUDY_FILE).write_text(_study_yaml(listed, equilibration))
    (root / "README.md").write_text(_readme(root.name))
    (root / "data.example.yaml").write_text(DATA_EXAMPLE)
    (root / ".gitignore").write_text(GITIGNORE)
    (root / "LICENSE-code").write_text(MIT.format(year=year, holder=holder))
    (root / "LICENSE-data").write_text(CC_BY.format(year=year, holder=holder))
    (root / "environment" / "README.md").write_text(_environment_readme())
    for folder, what in FOLDERS.items():
        readme = root / folder / "README.md"
        if folder != "environment" and not readme.exists():
            readme.write_text(f"# {folder}/\n\n{what}.\n")
    if git:
        created.commit = init_repository(root, "Create study with polyzymd study init")
    return created


def new_condition(folder: Path, runs: Path) -> None:
    """Write a template ``config.yaml`` and ``structures/`` into the new ``folder``, to fill in.

    The config writes its runs into ``runs``, given relative to the config.
    ``structures/`` holds placeholder files that say which inputs to add.
    """
    from polyzymd.core.branding import prepend_file_header
    from polyzymd.utils.templates import render_package_template

    (folder / "structures").mkdir(parents=True)
    config = render_package_template(
        "polyzymd.templates",
        "config_template.yaml",
        {
            "project_name": folder.name,
            "runs_directory": os.path.relpath(runs, folder.resolve()),
        },
    )
    (folder / "config.yaml").write_text(prepend_file_header(config, comment_prefix="#"))
    for name in ("protein", "ligand"):
        text = render_package_template("polyzymd.templates", f"{name}_placeholder.txt.jinja")
        (folder / "structures" / f"place_{name}_here.placeholder.txt").write_text(
            prepend_file_header(text, comment_prefix="#")
        )


def add_condition(
    study: str | Path,
    label: str,
    *,
    config: Path | None = None,
    new: bool = False,
    source: str | None = None,
    settings: Sequence[tuple[str, Any]] = (),
    factors: Mapping[str, Any] | None = None,
) -> Path:
    """Add a condition to an existing study folder and list it in ``study.yaml``.

    With ``config``, the config and its input files are copied as
    :func:`copy_condition` does; with ``source``, the config of that
    condition of the study is copied the same way, with ``name:`` set to the
    new condition's folder name, runs of its own and the source's
    ``factors``; with ``new``,
    ``conditions/<name>/`` gets a template config to fill in
    (:func:`new_condition`). ``settings``, pairs of a dotted config key and a
    value, are then set in a copied config (:func:`set_config_values`), and
    the result must load; ``factors`` set the condition's factors in
    ``study.yaml``, over those copied from ``source``. The new config writes its runs into
    :func:`runs_folder`. The condition is added as one line under
    ``conditions:``, so the rest of ``study.yaml``, comments included, is kept
    as written, and ``runs/`` is added to the ``.gitignore`` beside the
    ``runs/`` folder if it lacks it (:func:`ignore_runs`). Returns the new
    condition's config path.

    Raises
    ------
    ProtocolError
        If not exactly one of ``config``, ``source`` and ``new`` is given,
        the label or its folder is taken or reserved (:func:`check_label`), ``source`` is not a condition of
        the study, the config cannot be copied, ``settings`` are given with
        ``new``, a factor is not a single value, or the edited config does
        not load.
    """
    from polyzymd.analyses.study_file import _condition, find_study_file

    if [config is not None, new, source is not None].count(True) != 1:
        raise ProtocolError(
            "Give one of a config to copy, a condition to copy (from) or new.",
            hint="polyzymd study add-condition LABEL --config path/to/config.yaml, "
            "--from OTHER_LABEL, or --new.",
        )
    check_label(label, "condition")
    if settings and new:
        raise ProtocolError(
            "--set edits a copied config, and --new writes a template to fill in.",
            hint="Fill in the template by hand, or copy a config with --from or --config.",
        )
    if any(isinstance(value, (Mapping, list)) for value in (factors or {}).values()):
        raise ProtocolError("A factor must be a single value.", hint="--factor sbma_fraction=0.5")
    file = find_study_file(study)
    # Only the existing labels are read, so a new study with no condition yet,
    # or one still being filled in, takes its first conditions.
    conditions = _read_yaml(file).get("conditions") or {}
    root = file.parent
    folder = root / "conditions" / condition_folder(label)
    if label in conditions or folder.exists():
        raise ProtocolError(
            f"{file} already has a condition {label!r} or a folder {folder}.",
            hint="Choose another label.",
        )
    if source is not None:
        if source not in conditions:
            raise ProtocolError(
                f"{file} has no condition {source!r} to copy.",
                hint=f"Use one of: {', '.join(map(str, conditions)) or 'none yet'}.",
            )
        config, copied = _condition(source, conditions[source], file)
        factors = {**copied, **(factors or {})}
    else:
        factors = dict(factors or {})
    # A failure removes the condition folder and restores study.yaml, so the
    # study is left as it was.
    text = file.read_text()
    try:
        if config is not None:
            # A copy of another condition of the study is renamed as --new names
            # a template, and its runs go into a folder of its own named as its
            # runs/ folder is; a config from outside keeps its name and runs.
            runs = runs_folder(root, label)
            if source is None:
                copy_condition(Path(config), folder, runs)
            else:
                home = runs_home(root) / "runs"
                copy_condition(
                    Path(config),
                    folder,
                    runs,
                    name=condition_folder(label),
                    own_folder=runs.relative_to(home).as_posix(),
                    source_folder=runs_folder(root, source).relative_to(home).as_posix(),
                )
        else:
            new_condition(folder, runs_folder(root, label))
        if settings:
            _set_and_check(folder / "config.yaml", settings)
        entry = str((folder / "config.yaml").relative_to(root))
        _list_entry(
            file, "conditions", label, {"config": entry, "factors": factors} if factors else entry
        )
        if label not in (_read_yaml(file).get("conditions") or {}):
            raise ProtocolError(
                f"{file}: the condition {label!r} was not listed under conditions:.",
                hint="Add it under conditions: by hand.",
            )
    except BaseException:
        shutil.rmtree(folder, ignore_errors=True)
        file.write_text(text)
        raise
    if config is not None and source is None:
        record_data_location(root, label, Path(config))
    ignore_runs(runs_home(root))
    return folder / "config.yaml"


def _set_and_check(target: Path, settings: Sequence[tuple[str, Any]]) -> None:
    """Set ``settings`` in ``target`` and check that it loads, explaining what does not."""
    from pydantic import ValidationError

    from polyzymd.config.loader import config_error_lines, load_config

    set_config_values(target, settings)
    try:
        load_config(target)
    except (OSError, ValueError) as exc:
        lines = config_error_lines(exc, target) if isinstance(exc, ValidationError) else []
        problems = "; ".join(text for kind, text in lines if kind == "error") or str(exc)
        raise ProtocolError(
            f"The config with --set does not load: {problems}",
            hint="Correct the --set keys or values; the study is left unchanged.",
        ) from exc


def _read_yaml(file: Path) -> dict:
    """Read a YAML mapping, or raise ProtocolError saying the file is not one."""
    import yaml

    try:
        raw = yaml.safe_load(file.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ProtocolError(
            f"Cannot read {file}: {exc}", hint="Check that it is valid YAML."
        ) from exc
    if not isinstance(raw, dict):
        raise ProtocolError(f"{file} must be a YAML mapping.", hint="See the documented example.")
    return raw


def _list_entry(file: Path, block: str, label: str, path: str | dict) -> None:
    """Insert ``label: path`` as the last entry of the top-level ``block:`` of ``file``, keeping the rest.

    A mapping, such as ``{config: ..., factors: {...}}``, is written on one line in flow style.
    """
    import json

    import yaml

    lines = file.read_text().splitlines()
    # A JSON string is a double-quoted YAML string, so any label and path read back as written.
    value = (
        yaml.safe_dump(
            path, default_flow_style=True, sort_keys=False, allow_unicode=True, width=10**6
        ).strip()
        if isinstance(path, dict)
        else json.dumps(path, ensure_ascii=False)
    )
    entry = f"  {json.dumps(label, ensure_ascii=False)}: {value}"
    for index, line in enumerate(lines):
        match = re.match(rf"^{block}:\s*(\{{\s*\}})?\s*(#.*)?$", line)
        if match:
            if match.group(1):
                lines[index] = f"{block}:"
                lines.insert(index + 1, entry)
                break
            end = index + 1
            while end < len(lines) and (
                lines[end].startswith((" ", "\t")) or not lines[end].strip()
            ):
                end += 1
            while end > index + 1 and not lines[end - 1].strip():
                end -= 1
            lines.insert(end, entry)
            break
    else:
        raise ProtocolError(
            f"{file} has no top-level '{block}:' block to add to.",
            hint=f"Add {label} under {block}: by hand.",
        )
    file.write_text("\n".join(lines) + "\n")
